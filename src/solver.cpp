// Use of this source code is governed by an MIT-style license that can be found
// in the LICENSE file or at https://opensource.org/licenses/MIT. 
//
// Copyright (c) 2022 Maxence Reberol

#include "solver.h"

#include <ctime>

/* Gecode includes */
#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/minimodel.hh>

namespace hbl {
  using namespace Gecode;
  using std::vector;
  using std::pair;

  double Cpu() {
    std::clock_t c_start = std::clock();
    double t = double(c_start) / CLOCKS_PER_SEC;
    return t;
  }

  /* factor on integers to get a "floating-point energy" but using
   * only integers. This is because the mixed-integer solve of Gecode
   * is not as stable as the integer one. */
  constexpr int FACTOR = 100; 

  struct PolygonReduction {
    vector<pair<int,int>> var_fixed;
    vector<pair<int,int>> var_shift;
  };

  bool reduce_polygon(const PolygonReduction& poly, 
      vector<PolygonReduction>& sub_polys,
      const vector<int>& valence_min) {
    sub_polys.push_back(poly);

    /* Ending condition */
    if (poly.var_shift.size() == 3) {
      return true;
    }

    /* Recursive call */
    for (int i = 0; i < poly.var_shift.size(); ++i) {
      int var = poly.var_shift[i].first;
      int shift = poly.var_shift[i].second;
      if ((int) valence_min[var] + shift > 1) continue;
      int iprev = (poly.var_shift.size()+i-1)%poly.var_shift.size();
      int inext = (poly.var_shift.size()+i+1)%poly.var_shift.size();
      PolygonReduction sub;
      sub.var_fixed = poly.var_fixed;
      sub.var_fixed.push_back({var,1-shift});
      for (int j = 0; j < poly.var_shift.size(); ++j) if (j != i) {
        if (j == iprev) {
          sub.var_shift.push_back({poly.var_shift[j].first,poly.var_shift[j].second-1});
        } else if (j == inext) {
          sub.var_shift.push_back({poly.var_shift[j].first,poly.var_shift[j].second-1});
        } else {
          sub.var_shift.push_back({poly.var_shift[j].first,poly.var_shift[j].second});
        }
      }
      reduce_polygon(sub, sub_polys, valence_min);
    }

    return true;
  }

  class FindEdgeValences : public IntMinimizeSpace {
    protected:
      IntVarArray n;
      IntVar objective;
      const size_t N;
      const vector<vector<int>>& polygons;
      const vector<double>& x_i;
      const vector<std::array<int,2> >& n_min_max;
      const std::vector<int>& n_priority;

    public:
      FindEdgeValences(size_t N_, 
          const vector<vector<int>>& polygons_,
          const vector<double>& x_ideal_,
          const vector<std::array<int,2> >& n_min_max_,
          const std::vector<int>& n_prio_
          ) 
        : n(*this, N_), objective(*this,0,N_*VAL_MAX*VAL_MAX*FACTOR*FACTOR), N(N_), polygons(polygons_),
        x_i(x_ideal_), n_min_max(n_min_max_), n_priority(n_prio_)
    {
      /* Variables */
      vector<int> valence_min(N,1);
      for (size_t i = 0; i < N; ++i) {
        n[i] = IntVar(*this,n_min_max[i][0],n_min_max[i][1]);
      }

      /* Objective function to minimize */
      IntVarArray squares(*this,n.size());
      for (size_t i = 0; i < N; ++i) {
        int n_optimal = (int) std::round(double(FACTOR) * x_i[i]);
        squares[i] = expr(*this, (FACTOR * n[i] - n_optimal) * (FACTOR * n[i] - n_optimal));
      }
      rel(*this, objective == sum(squares));

      /* Loop over polygons and create integer constraints */
      size_t nbc = 0;
      size_t nbs = 0;
      for (size_t i = 0; i < polygons.size(); ++ i) {
        const vector<int>& vert = polygons[i];

        /* Loop over possible polygon reductions */
        PolygonReduction poly;
        for (size_t j = 0; j < vert.size(); ++j) {
          poly.var_shift.push_back({vert[j],0});
        }
        vector<PolygonReduction> sub_polys;
        reduce_polygon(poly, sub_polys, valence_min);
        for (size_t j = 0; j < sub_polys.size(); ++j) {
          const PolygonReduction& cpoly = sub_polys[j];
          int nfixed = cpoly.var_fixed.size();
          int nfree = cpoly.var_shift.size();

          /* Fixed variables during the reduction */
          IntVarArray fixed_vals(*this,nfixed);
          IntVar n_respected(*this,0,nfixed+1);
          if (nfixed > 0) {
            for (int k = 0; k < nfixed; ++k) {
              int v = (int) cpoly.var_fixed[k].first;
              int s = cpoly.var_fixed[k].second;
              fixed_vals[k] = expr(*this, n[v] - s);
            }
            count(*this,fixed_vals,0,IRT_EQ,n_respected);
          }

          /* Remaining variables in the sub poly */
          IntVarArray vals(*this,cpoly.var_shift.size());
          for (size_t k = 0; k < cpoly.var_shift.size(); ++k) {
            int v = (int) cpoly.var_shift[k].first;
            int s = cpoly.var_shift[k].second;
            vals[k] = expr(*this, n[v] + s);
          }
          IntVar n_val1(*this,0,nfree+1);
          IntVar n_val2(*this,0,nfree+1);
          IntVar n_valN(*this,0,nfree+1);
          count(*this,vals,1,IRT_EQ,n_val1);
          count(*this,vals,2,IRT_EQ,n_val2);
          count(*this,vals,0,IRT_LQ,n_valN);


          if (nfixed > 0) {
            /* Half reify */
            BoolVar cond1 = expr(*this, n_respected == nfixed);
            Reify r1(cond1, RM_IMP);
            if (nfree == 3) { /* triangle */
              rel(*this,n_val1,IRT_NQ,1,r1);
              rel(*this,n_val1,IRT_NQ,2,r1);
              nbc += 2;
            }
            rel(*this,n_val2,IRT_NQ,nfree-1,r1);
            rel(*this,n_valN,IRT_EQ,0,r1);
            nbc += 2;
          } else {
            if (nfree == 3) { /* triangle */
              rel(*this,n_val1,IRT_NQ,1);
              rel(*this,n_val1,IRT_NQ,2);
              nbc += 2;
            }
            rel(*this,n_val2,IRT_NQ,nfree-1);
            rel(*this,n_valN,IRT_EQ,0);
            nbc += 2;
          }
        }
      } /* end of loop on polygons */

      printf("Gecode space stats: %li integer variables, %li trgl. constraints (from %li polygons)\n", 
          N, nbc, polygons.size());

      /* Set the branching functions */
      branch(*this, n, INT_VAR_MERIT_MAX(&choose_variable), INT_VAL(&choose_value));
    }

      FindEdgeValences(FindEdgeValences& s) : IntMinimizeSpace(s), N(s.N), objective(s.objective), polygons(s.polygons), 
      x_i(s.x_i), n_min_max(s.n_min_max), n_priority(s.n_priority) {
        n.update(*this, s.n);
        objective.update(*this, s.objective);
      }

      virtual IntMinimizeSpace* copy(void) {
        return new FindEdgeValences(*this);
      }

      virtual IntVar cost(void) const {
        return objective;
      }

      double get_cost() const {
        return double(objective.val()) / double(FACTOR);
      }

      void print_cost(void) const {
        printf("cost: %.2f\n", double(objective.val()) / double(FACTOR));
      }

      void print(void) const {
        std::cout << n << std::endl;
      }

      void extract_solution(vector<int>& slt) {
        slt.resize(n.size());
        for (size_t i = 0; i < slt.size(); ++i) {
          slt[i] = (int) n[i].val();
        }
      }

      /* return closest value to x_ideal in input domain */
      static int choose_value(const Space& home, IntVar x, int i) {
        const FindEdgeValences& feb = static_cast<const FindEdgeValences&>(home);
        double ideal = feb.x_i[i];
        int val = -1;
        double dmin = DBL_MAX;
        IntVarValues values(x);
        while(values()) {
          double d = std::abs(double(values.val()) - ideal);
          if (d < dmin) {
            val = values.val();
            dmin = d;
          }
          ++values;
        }
        return val;
      }

      /* return priority for variable, using branching variable selection */
      static int choose_variable(const Space& home, IntVar x, int i) {
        const FindEdgeValences& feb = static_cast<const FindEdgeValences&>(home);
        return int(feb.n_priority[i]);
      }
  };

  class CustomStop : public Search::Stop {
    public:
      bool use_second;
    protected:
      Search::TimeStop* ts1;
      Search::TimeStop* ts2;
    public:
      /// Initialize stop object
      CustomStop(unsigned int t1, unsigned t2) 
        : ts1(new Search::TimeStop(t1)), ts2(new Search::TimeStop(t2)), use_second(false) {}
      /// Test whether search must be stopped
      virtual bool stop(const Search::Statistics& s, const Search::Options& o) {
        return use_second ? ts2->stop(s,o) : ts1->stop(s,o);
      }
      /// Destructor
      ~CustomStop(void) {
        delete ts1; 
        delete ts2;
      }
  };

  bool solveAllHexLayerTopology(
      size_t N,
      const std::vector<std::vector<int> >& polygons,
      const std::vector<double>& x_i,
      const std::vector<std::array<int,2> >& n_min_max,
      const std::vector<int>& n_priority,
      std::vector<int>& n,
      bool& stopped,
      std::vector<double>& iterTime,
      double timeMaxInit,
      double timeMaxImprove) 
  {
    FindEdgeValences m(N, polygons, x_i, n_min_max, n_priority);
    unsigned long int tl_init_ms = timeMaxInit;
    unsigned long int tl_2d_ms = timeMaxImprove;
    Search::Options opt;
    Search::TimeStop ts(tl_init_ms);
    CustomStop cstop(tl_init_ms, tl_2d_ms);
    opt.stop = &cstop;
    printf("Branch and bound search (max time init: %.1fs, improve: %.1fs) ... \n",
        double(tl_init_ms)/1.e3, double(tl_2d_ms)/1.e3);

    double t0 = Cpu();
    BAB<FindEdgeValences> e(&m,opt);
    size_t nf = 0;
    bool found = false;
    double last_cost = 0;
    while (FindEdgeValences* s = e.next()) {
      s->extract_solution(n);
      found = true;
      double cost = s->get_cost();
      double t = Cpu()-t0;
      printf("(%.1f, %.1fsec), ", cost, t);
      fflush(stdout);
      last_cost = cost;
      cstop.use_second = true;
      nf += 1;
      iterTime.push_back(t);
      delete s;
    }
    stopped = e.stopped();
    double t = Cpu()-t0;
    printf("  done (nIter: %li, stopped: %i, %.1fsec, found: %i)\n", nf, (int) stopped, t, int(found));

    return found;
  }

  bool solveAllHexLayerTopologyAllSolutions(
      size_t N,
      const std::vector<std::vector<int> >& polygons,
      const std::vector<double>& x_i,
      const std::vector<std::array<int,2> >& n_min_max,
      const std::vector<int>& n_priority,
      std::vector<std::vector<int> >& solutions) {

    FindEdgeValences m(N, polygons, x_i, n_min_max, n_priority);
    printf("DFS search to find all solutions ... \n");
    DFS<FindEdgeValences> e(&m);
    while (FindEdgeValences* s = e.next()) {
      vector<int> solution;
      s->extract_solution(solution);
      solutions.push_back(solution);
      delete s;
    }
    printf("found %li solutions \n", solutions.size());

    return solutions.size() > 0;
  }
}
