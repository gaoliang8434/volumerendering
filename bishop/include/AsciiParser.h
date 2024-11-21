//-----------------------------------------------------------
//
//  AsciiParser.h
//
//  Parses ascii text files, separates into tokens, and classifies
//  the tokens by type.
//
//  Copyright (c) 2005, Finelight Visual Technology, Inc.
//
//------------------------------------------------------------


#ifndef ____LUX_ASCIIPARSER_H____
#define ____LUX_ASCIIPARSER_H____

#include <string>
#include <vector>
using namespace std;

namespace lux
{

class AsciiParser
{
  public:

    AsciiParser();
   ~AsciiParser(){}

    //! Reads the file and parses contents into tokens.
    /*!
     * Returns true if parse successful, false otherwise
     */
    const bool ParseFile( const string& filename );

    //! Restart passing tokens from the beginning
    void Rewind();

    //! Advance to next token.  Return false if move past the last token.
    const bool GetToken();

    // make sure we are pointing to a valid token.
    const bool IsToken() const;

    //! Return true if the token is text
    const bool IsText()    const;
    //! Return true if the token is an integer
    const bool IsInteger() const;
    //! Return true if the token is a float (or double)
    const bool IsFloat()   const;
    //! Return true if the token is not recognized.
    const bool IsSeparator()   const;
    //! Return true if the token is a separator.
    const bool IsOther()   const;
    //! Return true if the token is an end of line.
    const bool IsEOL()   const;

    //! Retrieve value of an integer token. Return 0 if the token is not an integer.
    const int IntegerValue()  const;
    //! Retrieve the value of a float token. Return 0 if the token is not a float.
    const double FloatValue() const;
    //! Retrieve the value of a text token.  Return "" iif the token is not text.
    const string& TextValue() const;

    void AddWhiteSpaceCharacter( const string& ws );

    const vector<string>& Tokens() const { return tokens; }

  private:

    enum TokenType { OTHER=0, STRING, INTEGER, FLOAT, SEPARATOR, EOL };
    vector<string>    tokens;
    vector<TokenType> token_type;
    int current_token; 
    const string blank;
    string whitespace_characters;
    const string valid_integer_characters;
    const string valid_float_characters;
    const string valid_separator_characters;

    // turn a line of input text into tokens
    void Tokenize( const string& line );
    void AddToken( const string& line );
    void AddEOLToken();
    void AddOtherToken();
};

}

#endif
