

#ifndef ___LOGGER_H___
#define ___LOGGER_H___


#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <sys/types.h>
#include <unistd.h>

using namespace std;
namespace lux
{



class Logger
{
  public:

    Logger()
    {
       setData("generic");
       writeMessage("OPEN");
    }

    Logger( const string logfile, int argc, char** argv )
    {
       setData( logfile );
       string msg = "OPEN ";
       for(int i=0;i<argc;i++ )
       {
          msg += string(argv[i]) + " ";
       }
       writeMessage( msg );
    }


   ~Logger()
    {
       writeMessage("CLOSE");
    }


  private:

    template <typename T> 
    string tostr(const T& t) { std::ostringstream os; os<<t; return os.str(); }


    string pid;
    string user;
    string baseTag;
    string logFileName;

    void setData(const string log)
    {
       user = getenv("USER");
       logFileName = "/scratch/" + user + "/projects/logs/" + log + ".log";
       pid = tostr(getpid()) + " ";
       baseTag = pid + " ";
    }

    const string currentTag()
    {
       time_t seconds = time(NULL);
       string timestamp = ctime( &seconds );
       return (baseTag + timestamp.substr(0,timestamp.size()-1) + " ");
    }

    const bool writeMessage( const string msg )
    {
       string out = currentTag() + msg + "\n";
       ofstream logit;
       logit.open( logFileName.c_str(), ios_base::app );
       if( logit.good() )
       {
          logit << out;
          logit.close();
	  return true;
       }
       else
       {
          logit.open( logFileName.c_str() );
	  if( logit.good() )
          {
             logit << out;
             logit.close();
	     return true;
          }
          cerr << "Unable to write log message \"" << msg << "\" to file " << logFileName << endl;
       }
       return false;
    }

};
}




#endif
