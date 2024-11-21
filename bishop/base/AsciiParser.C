//-----------------------------------------------------------
//
//  AsciiParser.C
//
//  Parses ascii text files, separates into tokens, and classifies
//  the tokens by type.
//
//------------------------------------------------------------

#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;
#include "AsciiParser.h"

namespace lux
{

AsciiParser::AsciiParser() :
   current_token           (-1),
   blank                   (""),
   whitespace_characters   (" \r\t\n"),
   valid_integer_characters("-0123456789"),
   valid_float_characters  ("-0123456789.e"),
   valid_separator_characters ("/")
{}   


void AsciiParser::AddWhiteSpaceCharacter( const string& ws )
{
   whitespace_characters.append(ws);
}
	
const bool AsciiParser::ParseFile( const string& filename )
{
   // Open file, read characters, parse into tokens, assemble token list.

   ifstream file( filename.c_str() );
   if( !file ){ return false; }

   string token_candidate;
   int count = 0;
   while( !file.eof() )
   {
      // get next line
      string line;
      getline(file,line);
      // split into tokens
      Tokenize( line ); 
   }
   current_token = -1;
   return true;
}

void AsciiParser::Rewind() { current_token = -1; }
const bool AsciiParser::GetToken()
{
   ++current_token;
   if( current_token < (int)tokens.size() )
   { 
      return true; 
   }
   else 
   { 
      return false; 
   }
}

const bool AsciiParser::IsText() const 
{
   if( !IsToken() ) { return false; }
   return ( token_type[current_token] == STRING );
}

const bool AsciiParser::IsInteger() const
{
   if( !IsToken() ) { return false; }
   return ( token_type[current_token] == INTEGER );
}

const bool AsciiParser::IsFloat() const
{
   if( !IsToken() ) { return false; }
   return ( token_type[current_token] == FLOAT );
}

const bool AsciiParser::IsSeparator() const
{
   if( !IsToken() ) { return false; }
   return ( token_type[current_token] == SEPARATOR );
}

const bool AsciiParser::IsEOL() const
{
   if( !IsToken() ) { return false; }
   return ( token_type[current_token] == EOL );
}

const bool AsciiParser::IsOther() const
{
   if( !IsToken() ) { return false; }
   return ( token_type[current_token] == OTHER );
}

const int AsciiParser::IntegerValue() const
{
   if( !IsInteger() ){ return 0; }
   return atoi( tokens[current_token].c_str() );
}

const double AsciiParser::FloatValue() const
{
   if( !IsFloat() ){ return 0; }
   return atof( tokens[current_token].c_str() );
}

const string& AsciiParser::TextValue() const
{
   if( !IsText() ){ return blank; }
   return tokens[current_token];
}
    
const bool AsciiParser::IsToken() const
{
   return ( ( current_token >= 0 ) && ( current_token < (int)tokens.size() ) );
}

void AsciiParser::Tokenize( const string& line )
{
   string tline(line);

   // remove leading and trailing whitespace
   string::size_type  notwhite = tline.find_first_not_of(whitespace_characters);
   tline.erase(0,notwhite);
   notwhite = tline.find_last_not_of(whitespace_characters);
   tline.erase(notwhite+1);

   // split up by separators
   vector<string> tlineset;
   string::size_type notseparator = tline.find_first_of(valid_separator_characters);
   while( notseparator != string::npos )
   {
      string token;
      token.append( tline, 0,notseparator );
      tlineset.push_back(token);
      tlineset.push_back( valid_separator_characters );
      tline.erase(0,notseparator+1);
      notseparator = tline.find_first_of(valid_separator_characters);
   }
   if(tline.size()>0){ tlineset.push_back( tline ); }
  

   // now find whitespace characters that separate tokens
   for( size_t ic=0;ic<tlineset.size();ic++)
   {
      tline = tlineset[ic];
      notwhite = tline.find_first_of(whitespace_characters);
      while( notwhite != string::npos )
      {
         string token;
         token.append( tline, 0,notwhite );
         AddToken( token );
         tline.erase(0,notwhite);
         // remove leading whitespace again
         notwhite = tline.find_first_not_of(whitespace_characters);
         tline.erase(0,notwhite);
         // find next white space
         notwhite = tline.find_first_of(whitespace_characters);
      }
      if(tline.size()>0){ AddToken( tline ); }
   }
   // add and eol token
   AddEOLToken();
}

void AsciiParser::AddToken( const string& token )
{
   tokens.push_back( token );
   // need to also figure out the kind of token
   if( token.find_first_not_of( valid_integer_characters ) == string::npos )
   {
      token_type.push_back( INTEGER );
   }
   else if( token.find_first_not_of( valid_float_characters ) == string::npos )
   { 
      token_type.push_back( FLOAT );
   }
   else if( token.find_first_not_of( valid_separator_characters ) == string::npos )
   { 
      token_type.push_back( SEPARATOR );
   }
   else
   {
      token_type.push_back( STRING );
   }
}

void AsciiParser::AddEOLToken()
{
   tokens.push_back( blank );
   token_type.push_back( EOL );
}

void AsciiParser::AddOtherToken()
{
   tokens.push_back( blank );
   token_type.push_back( OTHER );
}

}

