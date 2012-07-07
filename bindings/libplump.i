%module libplump
%feature("nestedworkaround");
%include "std_vector.i"
%include "std_string.i"

namespace std {
   %template(VectorInt) vector<int>;
   %template(VectorDouble) vector<double>;
}

%{
/* Includes the header in the wrapper code */
#include "libplump/hpyp_model.h"
#include "libplump/config.h"
#include "libplump/utils.h"
#include "libplump/node_manager_interface.h"
#include "libplump/node_manager.h"
#include "libplump/context_tree.h"
#include "libplump/hpyp_restaurant_interface.h"
#include "libplump/hpyp_restaurants.h"
#include "libplump/hpyp_parameters_interface.h"
#include "libplump/hpyp_parameters.h"
#include "libplump/random.h"
#include "libplump/serialization.h"
#include "libplump/pyp_sample.h"
#include "libplump/stirling.h"
%}

namespace gatsby { namespace libplump {

class DFSPathIterator {
public:
  DFSPathIterator(gatsby::libplump::ContextTree::NodeId root, const gatsby::libplump::INodeManager& nm, 
                  const gatsby::libplump::ContextTree& ct);
private:
  DFSPathIterator();
};

typedef int int32_t;

%nestedworkaround gatsby::libplump::ContextTree::DFSPathIterator;



}}

%ignore getDFSPathIterator;

/* Parse the header file to generate wrappers */
%include "libplump/config.h"
%include "libplump/utils.h"
%include "libplump/serialization.h"
%include "libplump/context_tree.h"
%include "libplump/node_manager_interface.h"
%include "libplump/node_manager.h"
%include "libplump/hpyp_restaurant_interface.h"
%include "libplump/hpyp_parameters_interface.h"
%include "libplump/hpyp_restaurants.h"
%include "libplump/hpyp_parameters.h"
%include "libplump/random.h"
%include "libplump/hpyp_model.h"
%include "libplump/pyp_sample.h"
%include "libplump/stirling.h"
 
%{
namespace gatsby { namespace libplump {
  typedef ContextTree::DFSPathIterator DFSPathIterator;
}}
%}

namespace gatsby { namespace libplump {
%template(pushCharFileToVec) pushFileToVec<unsigned char, seq_type>;
%template(pushIntFileToVec) pushFileToVec<int, seq_type>;
%template(prob2loss) prob2loss<double>;
}}

%init %{
  gatsby::libplump::init_rng();
%}


