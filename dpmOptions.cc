#define WANT_STREAM
#define WANT_MATH

#include "dpmOptions.h"

using namespace Utilities;

namespace DPM {

dpmOptions* dpmOptions::gopt = NULL;

  void dpmOptions::parse_command_line(int argc, char** argv, Log& logger)
  {
    // do once to establish log directory name
    for(int a = options.parse_command_line(argc, argv); a < argc; a++);
    
    
    if(help.value() || ! options.check_compulsory_arguments())
      {
	options.usage();
	//throw Exception("Not all of the compulsory arguments have been provided");
	exit(2);
      }

    
  }
  
}
