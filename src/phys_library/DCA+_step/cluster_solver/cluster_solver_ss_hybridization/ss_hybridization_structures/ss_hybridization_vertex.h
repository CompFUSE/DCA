//-*-C++-*-

#ifndef SS_HYBRIDIZATION_VERTEX_H
#define SS_HYBRIDIZATION_VERTEX_H

namespace DCA 
{  
  /*! 
   *  \brief   This class organizes a vertex, stores the time points of a \f$c\f$ and \f$c^{\dagger}\f$.
   *  \author  Peter Staar
   *  \author  Bart Ydens
   *  \version 1.0
   */
  class Hybridization_vertex
  {
    typedef Hybridization_vertex this_type;

  public:

    Hybridization_vertex();
    Hybridization_vertex(double t_start, double t_end);

    double t_start() const;
    double t_end()   const;

    void set_t_start(double t_start);
    void set_t_end  (double t_end);

    this_type& operator=(      this_type& other_vertex);
    this_type& operator=(const this_type& other_vertex);
    
  private:

    double t_start_val, t_end_val;
  };

  Hybridization_vertex::Hybridization_vertex():
    t_start_val(0),
    t_end_val(0)
  {}

  Hybridization_vertex::Hybridization_vertex(double t_st, double t_e):
    t_start_val(t_st),
    t_end_val(t_e)
    {}

  double Hybridization_vertex::t_start() const 
  {
    return t_start_val;
  }
  
  double Hybridization_vertex::t_end() const 
  {
    return t_end_val;
  }

  void Hybridization_vertex::set_t_start(double t_st)
  {
    t_start_val = t_st;
  }

  void Hybridization_vertex::set_t_end(double t_e) 
  {
    t_end_val = t_e;
  }


  bool operator<(const Hybridization_vertex& t1, const Hybridization_vertex& t2) {
    return t1.t_start() < t2.t_start();
  }
  
  bool operator<(const Hybridization_vertex& t1, const double t2) {
    return t1.t_start() < t2;
  }
  
  bool operator>(Hybridization_vertex t1, Hybridization_vertex t2) {
    return t1.t_start() > t2.t_start();
  }
  
  bool operator==(Hybridization_vertex t1, Hybridization_vertex t2) {
    return t1.t_start() == t2.t_start();
  }

  Hybridization_vertex& Hybridization_vertex::operator=(Hybridization_vertex& other_vertex)
  {
    t_start_val = other_vertex.t_start();
    t_end_val   = other_vertex.t_end();
    
    return *this;
  }

  Hybridization_vertex& Hybridization_vertex::operator=(const Hybridization_vertex& other_vertex)
  {
    t_start_val = other_vertex.t_start();
    t_end_val   = other_vertex.t_end();

    return *this;
  }

}

#endif

