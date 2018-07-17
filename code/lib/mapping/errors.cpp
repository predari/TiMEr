#include "errors.h"

/******************************************************/
/*                        error                       */
/******************************************************/
void error(std::string errorMessage, int errorCode) {
  std::cout << "ERROR: " << errorMessage << "\n";
  exit(errorCode);
}
