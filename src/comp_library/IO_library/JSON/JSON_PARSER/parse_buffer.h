//-*-C++-*-

#ifndef  JSONPARSER_PARSEBUFFER_H
#define  JSONPARSER_PARSEBUFFER_H

namespace IO
{
  namespace JSONPARSER
  {
    class ParseBuffer// : 
    //     public type_mixing
    {
    public:

      ParseBuffer();
      ~ParseBuffer();

      void clear();
    
      void put(wchar_t wc);

      std::wstring str(); 

      std::string to_string(); 

    public:

      std::vector<wchar_t> theCharacters;
      bool                 trace;    
    };

    ParseBuffer::ParseBuffer():
      //     type_mixing(),
      theCharacters(),
      trace(false)
    {}

    ParseBuffer::~ParseBuffer()
    {}

    void ParseBuffer::clear() 
    {
      theCharacters.clear();
    
      if(trace) 
	std::wcout << L"   ParseBuffer: clear()\n";
    }
    
    void ParseBuffer::put(wchar_t wc) 
    {
      theCharacters.push_back(wc);
    }

    std::wstring ParseBuffer::str() 
    {
      return std::wstring(theCharacters.begin(), theCharacters.end());
    }

    std::string ParseBuffer::to_string() 
    {
      return std::string(theCharacters.begin(), theCharacters.end());
    }

  }

}

#endif
