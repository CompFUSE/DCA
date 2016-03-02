//-*-C++-*-

#ifndef JSONPARSER_WHATEVER_H
#define JSONPARSER_WHATEVER_H

namespace IO
{
  namespace JSONPARSER
  {
    class Whatever 
    {
    public:
    
      typedef std::map<std::wstring,Whatever> WhateverMap;
      typedef std::vector<Whatever>           WhateverVector;

    public:

      Whatever();
      ~Whatever();

      Whatever(JSON_whatever_type t);

      size_t size() const; 

      //std::wstring name()  const ;
      std::string name() const ;

      Whatever& operator[] (const std::wstring key);
      Whatever& operator[] (const std::string  key);
      Whatever& operator[] (size_t             index); 
    
      const Whatever& operator[] (const std::wstring key)   const ;
      const Whatever& operator[] (const std::string  key)   const ;
      const Whatever& operator[] (size_t             index) const ;

      static Whatever null();
      Whatever& back(); 
    
      template<typename T>
      Whatever& push_back(T& value); 
    
      Whatever& push_back();     
      Whatever& push_back_null();

      template<typename T>
      Whatever& operator = (const T& value);

    private:
    
      static std::wstring typeName(JSON_whatever_type t);
  
      static std::string ntypeName(JSON_whatever_type t);

      void collectName(std::wostringstream& nam) const;
    
      int count(std::wstring key); 
    
      template<typename T>
      void get(T& value); 
    
      void setNull(); 

      void set(const std::wstring& value); 
    
      void set(bool value); 
    
      void set(int value); 
    
      void set(double value); 

      //private:
    public:

      JSON_whatever_type type;
      Whatever*          parent;

      std::wstring   valueString;

      WhateverMap    whateverMap;
      WhateverVector whateverVector;
      bool           whateverBool;
      int            whateverInteger;
      double         whateverDouble;

      std::wstring   myKey;
      int            myIndex;

      std::string    filename;
      int            startPos;
      int            endPos;
    };
  
    Whatever::Whatever():
      type(WHATEVER_UNKNOWN),
      parent(0),
      valueString(),
      whateverMap(),
      whateverVector(),
      whateverBool(true),
      whateverInteger(0),
      whateverDouble(0),
      myKey(L"?"),
      myIndex(-1),
      startPos(0),
      endPos(0)
    {}
  
    Whatever::Whatever(JSON_whatever_type t):
      type(t),
      parent(0),
      valueString(),
      whateverMap(),
      whateverVector(),
      whateverBool(true),
      whateverInteger(0),
      whateverDouble(0),
      myKey(L"?"),
      myIndex(-1),
      startPos(0),
      endPos(0)
    {}

    Whatever::~Whatever()
    {}
  
    std::wstring Whatever::typeName(JSON_whatever_type t) 
    {
      switch (t) 
	{
	case WHATEVER_MAT:
	  return L"WHATEVER_MAT";
	case WHATEVER_MAP:
	  return L"WHATEVER_MAP";
	case WHATEVER_VECTOR:
	  return L"WHATEVER_VECTOR";
	case WHATEVER_MATRIX:
	  return L"WHATEVER_MATRIX";
	case WHATEVER_STRING:
	  return L"WHATEVER_STRING";
	case WHATEVER_INTEGER:
	  return L"WHATEVER_INTEGER";
	case WHATEVER_DOUBLE:
	  return L"WHATEVER_DOUBLE";
	case WHATEVER_BOOL:
	  return L"WHATEVER_BOOL";
	case WHATEVER_UNKNOWN:
	  return L"WHATEVER_UNKNOWN";
	default:
	  throw std::logic_error("Whatever::typeName given wrong type");
	}
    }
  
    std::string Whatever::ntypeName(JSON_whatever_type t) 
    {
      std::wstring wname(typeName(t));
      return std::string(wname.begin(), wname.end());
    }

    Whatever Whatever::null() 
    {
      Whatever result;
      result.setNull();
      return result;
    }
  
    std::string Whatever::name() const 
    {
      std::wostringstream nam;

      collectName(nam);
      nam << "{" << typeName(type) << "}";

      std::wstring wname = nam.str();//(name());

      return std::string(wname.begin(),wname.end());
    }

    void Whatever::collectName(std::wostringstream& nam) const 
    {
      if(parent == 0) 
	{
	  nam << "root";
	  return;
	}

      parent->collectName(nam);

      switch (parent->type) 
	{
	case WHATEVER_MAP:
	  nam << L"[" << myKey << L"]";
	  break;
	case WHATEVER_VECTOR:
	  nam << L"[" << myIndex << L"]";
	  break;
	default:
	  nam << L"[ ?? not vector or map! ]";
	}
    }
    
  
    Whatever& Whatever::operator[] (const std::string key) 
    {
      std::wstring wKey(key.begin(),key.end());

      assert(type == WHATEVER_MAP and whateverMap.count(wKey) == 1);

      return whateverMap[wKey];
    }
    
    Whatever& Whatever::operator[] (const std::wstring key) 
    {
      assert(type == WHATEVER_MAP and whateverMap.count(key) == 1);

      return whateverMap[key];
    }
  

    Whatever& Whatever::operator[] (size_t index) 
    {
      assert(type == WHATEVER_VECTOR);
      assert(index<whateverVector.size());

      return whateverVector[index];
    }


    const Whatever& Whatever::operator[] (const std::wstring key) const 
    {
      assert(type == WHATEVER_MAP and whateverMap.count(key) == 1);

      WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);

      return wm[key];
    }
        
    const Whatever& Whatever::operator[] (const std::string key) const 
    {
      std::wstring wKey(key.begin(),key.end());

      WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);

      return wm[wKey];
    }


    const Whatever& Whatever::operator[] (size_t index) const 
    {
      assert(type == WHATEVER_VECTOR);
      assert(index<whateverVector.size());
    
      return whateverVector[index];
    }
  
    size_t Whatever::size() const 
    {
      assert(type == WHATEVER_VECTOR);
      return whateverVector.size();
    }
  
    int Whatever::count(std::wstring key) 
    {    
      assert(type == WHATEVER_VECTOR);
      return whateverMap.count(key);
    }
  
    Whatever& Whatever::back() 
    {    
      assert(type == WHATEVER_VECTOR);
      return whateverVector.back();
    }
  
    template<typename T>
    Whatever& Whatever::push_back(T& value) 
    {    
      assert(type == WHATEVER_VECTOR);

      whateverVector.push_back(Whatever());
      whateverVector.back() = value;

      return *this;
    }
  
    Whatever& Whatever::push_back() 
    {
      assert(type == WHATEVER_VECTOR);
    
      whateverVector.push_back(Whatever());

      return *this;
    }
  
    Whatever& Whatever::push_back_null()
    {
      assert(type == WHATEVER_VECTOR);
    
      whateverVector.push_back(Whatever());
      whateverVector.back().setNull();

      return *this;
    }
  
    template<typename T>
    void Whatever::get(T& value) 
    {
      value = valueString;
    }
  
    void Whatever::setNull() 
    {
      type = WHATEVER_NULL;
      valueString = L"NULL";
    }
    
    template<typename T>
    Whatever& Whatever::operator=(const T& value) 
    {
      set(value);
      return *this;
    }

    void Whatever::set(const std::wstring& value) 
    {
      type = WHATEVER_STRING;
      valueString = value;
    }
  
    void Whatever::set(bool value) 
    {
      type = WHATEVER_BOOL;
      whateverBool = value;
    }
  
    void Whatever::set(int value) 
    {
      type = WHATEVER_INTEGER;
      whateverInteger = value;
    }
  
    void Whatever::set(double value) 
    {
      type = WHATEVER_DOUBLE;
      whateverDouble = value;
    }
    
  } 

}

#endif

