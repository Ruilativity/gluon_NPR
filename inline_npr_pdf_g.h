// -*- C++ -*-
/*! \file
 * \brief Inline construction of Landau gauge propagator
 *
 * Landau gauge Propagators for NPR  calculations
 */

#ifndef __inline_npr_h__
#define __inline_npr_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "handle.h"
#include "state.h"

using namespace QDP;

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGluonNprPDFEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineGluonNprPDFParams
  {
    InlineGluonNprPDFParams();
    InlineGluonNprPDFParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    int max_mom2 ; // max p^2
	int lmax ;
//    std::string output_type ;
    std::string filename ;

    struct NamedObject_t
    {
      std::string     gauge_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineGluonNprPDF : public AbsInlineMeasurement
  {
  public:
    ~InlineGluonNprPDF() {}
    InlineGluonNprPDF(const InlineGluonNprPDFParams& p) : params(p) {}
    InlineGluonNprPDF(const InlineGluonNprPDF& p) : params(p.params) {}
 
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 
    
  private:
    InlineGluonNprPDFParams params;
  };

}

#endif
