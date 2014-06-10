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
	  {
	    lhs = static_cast<bool>(w.whateverInteger);
	    return lhs;
	  }

	case WHATEVER_BOOL:
	  {
	    lhs = w.whateverBool;
	    return lhs;
	  }

	default: 
	  {
	    //std::cout << "double d <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a bool!\n";
	    throw std::logic_error(__FUNCTION__);
	  }
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
	    //std::cout << "int d <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a int!\n";
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
	    //std::cout << "float d <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a float!\n";
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
	    //std::cout << "double d <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a double!\n";
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
	    //std::cout << "std::complex<T>  <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a std::complex<T>!\n";
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
	    ////std::cout << "std::string d <= " << w.name() << " produced a type error!\n";
	    ////std::cout << " trying to assign a " << name(w.type) << " to a std::string!\n";
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
	    //std::ostringstream msg;
	    //msg << "std::vector<t> <= " << w.name() << " produced a type error!\n";
	    //msg << " trying to assign a " << name(w.type) << " to a std::vector<t>!\n";
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
	    //std::cout << "std::vector<t> <= " << w.name() << " produced a type error!\n";
	    //std::cout << " trying to assign a " << name(w.type) << " to a std::vector<t>!\n";
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

    /*
      template<typename StreamType>
      StreamType& operator << (StreamType& os, const Whatever& w) 
      {
      static int level = -1;

      typedef Whatever::WhateverMap        WhateverMap;
      typedef WhateverMap::const_iterator  WhateverMapItr;
      typedef std::string                  StringType;
    
      switch (w.type) 
      {
      case WHATEVER_MAT: 
      {
      std::string wfilename(w.filename.begin(),w.filename.end());
      os << "{ 'fileName': '" << wfilename << "'"
      << ", 'startPos': "  << w.startPos 
      << ", 'endPos': "    << w.endPos << "}";
      break;
      }

      case WHATEVER_MAP:
      {
      level += 1;

      os << "\n";
      for(int l=0; l<level; l++)
      os << "\t";
      os << "{\n";

      int index = 0;
      for(WhateverMapItr itr = w.whateverMap.begin(); itr != w.whateverMap.end(); itr++) 
      {
      const std::wstring& wkey = (*itr).first;
      const std::string    key(wkey.begin(), wkey.end());

      for(int l=0; l<level; l++)
      os << "\t";

      os << "\"" << key << "\" : " << (*itr).second; 

      if(w.whateverMap.size() == index+1)
      os << "";
      else
      os << ",\n";

      index += 1;
      }

      os << "\n";
      for(int l=0; l<level; l++)
      os << "\t";
      os << "}";


      level -= 1;

      break;
      }
	
      case WHATEVER_VECTOR:
      {
      os << "[";
      for(size_t i=0; i<w.whateverVector.size(); i++)
      {
      os << w.whateverVector[i];
      if(i<w.whateverVector.size()-1)
      os << ", ";
      }

      os << "]";
	  
      break;
      }

      case WHATEVER_MATRIX:
      os << "WHATEVER_MATRIX";
      break;

      case WHATEVER_STRING:
      {
      const std::string tmp(w.valueString.begin(), w.valueString.end());
      os << "\"" << tmp << "\"";
      }
      break;

      case WHATEVER_INTEGER:
      os << w.whateverInteger;
      break;

      case WHATEVER_DOUBLE:
      os << w.whateverDouble;
      break;

      case WHATEVER_BOOL:
      os << w.whateverBool;
      break;

      case WHATEVER_UNKNOWN:
      os <<"WHATEVER_UNKNOWN";
      break;

      default:
      throw std::logic_error("typeName given wrong type");
      }

      return os;
      }
    */
  
  }

}

#endif
