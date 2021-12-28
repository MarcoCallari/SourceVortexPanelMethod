#pragma once
#include <gsl/gsl_integration.h>

/*typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;*/
namespace gsl {
  template <size_t N>
  class StackWorkspace {
  public:
    StackWorkspace();
    gsl_integration_workspace* ptr();
  private:
    gsl_integration_workspace m_w;
    double m_alist[N];
    double m_blist[N];
    double m_rlist[N];
    double m_elist[N];
    size_t m_order[N];
    size_t m_level[N];
  };

  template<size_t N>
  inline StackWorkspace<N>::StackWorkspace() {
    m_w.limit = N;
    m_w.size = 0;
    m_w.nrmax = 0;
    m_w.i = 0;
    m_w.maximum_level = 0;
    m_w.alist = m_alist;
    m_w.blist = m_blist;
    m_w.rlist = m_rlist;
    m_w.elist = m_elist;
    m_w.order = m_order;
    m_w.level = m_level;
  }

  template<size_t N>
  inline gsl_integration_workspace* StackWorkspace<N>::ptr() {
    return &m_w;
  }

}