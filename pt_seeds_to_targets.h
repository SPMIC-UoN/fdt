/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <fstream>
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "meshclass/meshclass.h"
#include "probtrackOptions.h"
#include "particle.h"
#include "tractvols.h"




void read_masks(vector<string>& masks, const string& filename);
void seeds_to_targets();
