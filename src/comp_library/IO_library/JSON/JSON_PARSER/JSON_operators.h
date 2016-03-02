//-*-C++-*-

#ifndef JSONPARSER_OPERATORS_H
#define JSONPARSER_OPERATORS_H

namespace IO
{
  namespace JSONPARSER
  {
    bool& operator <= (bool& lhs, const Whatever& w) 
    {
      switch(w.type) 
	{	
	case WHATEVER_INTEGER:
	  
	  lhs = static_cast<bool>(w.whateverInteger);
	  return lhs;
	  
	case WHATEVER_BOOL:
      
	  lhs = w.whateverBool;
	  return lhs;
      
	default:    throw std::logic_error(__FUNCTION__);
	}
    }

    int& operator <= (int& lhs, const Whatever& w) 
    {
      switch(w.type) 
	{
	case WHATEVER_INTEGER: 
	  {
	    lhs = w.whateverInteger;
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }
  
    float& operator <= (float& lhs, const Whatever& w) 
    {
      switch (w.type) 
	{
	case WHATEVER_INTEGER:
	  {
	    lhs = static_cast<float>(w.whateverInteger);
	    return lhs;
	  }

	case WHATEVER_DOUBLE:
	  {
	    lhs = static_cast<float>(w.whateverDouble);
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }

    double& operator <= (double& lhs, const Whatever& w) 
    {
      switch (w.type) 
	{
	case WHATEVER_INTEGER:
	  {
	    lhs = static_cast<double>(w.whateverInteger);
	    return lhs;
	  }
	case WHATEVER_DOUBLE:
	  {
	    lhs = w.whateverDouble;
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }

    template<typename T>
    std::complex<T>& operator <= (std::complex<T>& lhs, const Whatever& w) 
    {
      switch(w.type) 
	{
	case WHATEVER_VECTOR: 
	  {
	    T v0;
	    v0 <= w.whateverVector[0];
	    T v1;
	    v1 <= w.whateverVector[1];
	    std::complex<T> result(v0,v1);
	    lhs = result;
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }

    std::string& operator <= (std::string& lhs, const Whatever& w) 
    {
      switch (w.type) 
	{
	case WHATEVER_STRING:
	  {
	    lhs = std::string(w.valueString.begin(), w.valueString.end());
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }

    template<typename T>
    std::vector<T>& operator <= (std::vector<T>& lhs, const Whatever& w) 
    {
      switch (w.type) 
	{
	case WHATEVER_VECTOR:
	  {
	    lhs.resize(w.whateverVector.size());
	    for (size_t i=0; i<lhs.size(); i++)
	      lhs[i] <= w.whateverVector[i];
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }

    template<typename T>
    std::vector<std::vector<T> >& operator <= (std::vector<std::vector<T> >& lhs, const Whatever& w) 
    {
      switch(w.type) 
	{
	case WHATEVER_VECTOR:
	  {
	    lhs.resize(w.whateverVector.size(),std::vector<T>(w[0].size()));
	    for (size_t i=0; i< lhs.size(); i++)
	      lhs[i] <= w.whateverVector[i];
	    return lhs;
	  }

	default: 
	  {
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }

    template<typename T>
    std::map<std::string,T>& operator <= (std::map<std::string,T>& lhs, const Whatever& w) 
    {
      switch (w.type) 
	{
	case WHATEVER_MAP: 
	  {
	    typedef std::map<std::wstring,Whatever>          WhateverMapType;
	    typedef typename WhateverMapType::const_iterator itr_type;

	    lhs.clear();

	    for(itr_type itr = w.whateverMap.begin(); itr != w.whateverMap.end(); itr++) {
	      const std::pair<std::wstring,Whatever>& pair(*itr);
	      std::string key(pair.first.begin(), pair.first.end());
	      lhs[key] <= pair.second;
	    }

	    return lhs;
	  }
	
	default: 
	  {
	    //std::cout << "std::map<t,t2> <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a std::map!\n";
	    throw std::logic_error(__FUNCTION__);
	  }
	}
    }
  
    template<typename T>
    ParseBuffer& operator >> (ParseBuffer& buffer, T& value) 
    {
      std::wistringstream is(buffer.str());
      is >> value;
      return buffer;
    }

  
  }//namespace JSONPARSER

}//namespace IO

#endif
