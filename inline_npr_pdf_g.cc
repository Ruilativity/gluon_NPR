/*! \file Id: inline_npr_pdf_w.cc  Modified by R. Zhang in Nov. 1st, 2019
 * \brief Inline construction of NPR propagator
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "io_general_class.h"
#include "inline_npr_pdf_g.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
  namespace InlineGluonNprPDFEnv
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
					      const std::string& path)
      {
	return new InlineGluonNprPDF(InlineGluonNprPDFParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "GLUON_NPR_PDF";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  } // end namespace


  //! NprPDF input
  void read(XMLReader& xml, const std::string& path, InlineGluonNprPDFParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! NprPDF output
  void write(XMLWriter& xml, const std::string& path, const InlineGluonNprPDFParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    //write(xml, "source_id", input.source_id);

    pop(xml);
  }


  // Param stuff
  InlineGluonNprPDFParams::InlineGluonNprPDFParams() { frequency = 0; }

  InlineGluonNprPDFParams::InlineGluonNprPDFParams(XMLReader& xml_in, const std::string& path)
  {
    try
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

       // Read in the output npr/source configuration info
      read(paramtop, "max_mom2", max_mom2);

      read(paramtop, "filename", filename);

      // Read in the output npr/source configuration info
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0)
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e)
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void
  InlineGluonNprPDFParams::write(XMLWriter& xml_out, const std::string& path)
  {
    push(xml_out, path);
    
    QDP::write(xml_out, "filename", filename);
    QDP::write(xml_out, "max_mom2", max_mom2);
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }



  // Function call
  void
  InlineGluonNprPDF::operator()(unsigned long update_no,
			       XMLWriter& xml_out)
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "npr");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }


  // Real work done here
  void
  InlineGluonNprPDF::func(unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();
    
    StopWatch snoop;
    snoop.reset();
    snoop.start();
    
    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
    catch( std::bad_cast )
      {
	QDPIO::cerr << InlineGluonNprPDFEnv::name << ": caught dynamic cast error"
		    << std::endl;
	QDP_abort(1);
      }
    catch (const std::string& e)
      {
	QDPIO::cerr << InlineGluonNprPDFEnv::name << ": std::map call failed: " << e
		    << std::endl;
	QDP_abort(1);
      }
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	  
	  
	  
    push(xml_out, "npr");
    write(xml_out, "update_no", update_no);
    XMLBufferWriter file_xml;
    push(file_xml, "GLUON_NPR_PDF");
    push(file_xml, "Output_version");
    write(file_xml, "out_version", 1);
    pop(file_xml);
    proginfo(file_xml);    // Print out basic program info
    params.write(file_xml, "Input");

    QDPIO::cout << InlineGluonNprPDFEnv::name << ": gluon npr calculation" << std::endl;
    
    proginfo(xml_out);    // Print out basic program info
    
    // Write out the input
    params.write(xml_out, "Input");
    
    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);
    
    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);
    write(xml_out, "Config_info", gauge_xml);

    write(file_xml, "Config_info", gauge_xml);
    pop(file_xml) ; //NPR_W

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);
    
    // Make sure that the source location is irrelevant
    SftMom phases(params.max_mom2, false,  -1); // 4D fourier transform
	std::vector<int> mom_serial(phases.numMom());

	  double g0=1.0;
	  LatticeColorMatrix trUmat;
	  multi1d<LatticeColorMatrix> ai(Nd);
	  for(int mu=0; mu<Nd; ++mu){
		  trUmat=2.0/Nc*imag(trace(u[mu]));
		  ai[mu]=1/(2*g0)*(u[mu]-adj(u[mu])-(trUmat-adj(trUmat))/Nc)/cmplx(Real(0.0),Real(1.0));
	  }
	  
	 multi1d<int> checkpoint(4);
	  for(int mu=0; mu<Nd; ++mu) checkpoint[mu]=Layout::lattSize()[mu]/2;
	  QDPIO::cout << "check A(x):("<< checkpoint[0] <<","<< checkpoint[1] <<","<< checkpoint[2] <<","<< checkpoint[3] <<")" <<  std::endl;
	  for(int mu=0; mu<Nd; ++mu){
		  QDPIO::cout << "A_"<< mu <<"(x):";
		  for(int i=0; i<Nc; ++i)
			  for(int j=0; j<Nc; ++j)
				  QDPIO::cout<< peekSite(ai[mu],checkpoint).elem().elem().elem(i,j)<< "\t";
		  QDPIO::cout<< std::endl;
	  } 
	  
	  

	multi2d<ColorMatrix> Ap(phases.numMom(),Nd);
	multi3d<Complex> Dp(phases.numMom(),Nd ,Nd);
	Complex shift_phase;
	for (int m=0; m < phases.numMom(); m++){
QDPIO::cout << "p: "<< phases.numToMom(m)[0]<<phases.numToMom(m)[1]<<phases.numToMom(m)[2]<<phases.numToMom(m)[3] << std::endl;
		mom_serial[m]=(50+phases.numToMom(m)[0])+(50+phases.numToMom(m)[1])*100+(50+phases.numToMom(m)[2]) *10000+(50+phases.numToMom(m)[3])*1000000;
		
		for(int mu=0; mu<Nd; ++mu){
			const Real twopi = 6.283185307179586476925286;
			Real p_dot_x;
			p_dot_x=phases.numToMom(m)[mu]*twopi/Layout::lattSize()[mu]/2.0;
			shift_phase=cmplx(cos(-p_dot_x),sin(-p_dot_x));
			Ap[m][mu] = shift_phase*sum(phases[m]*ai[mu]);
		}
	}
	for (int m=0; m < (1+phases.numMom())/2; m++){
		for(int mu=0; mu<Nd; ++mu)
			for(int nu=0; mu<Nd; ++nu){
				Dp[m][mu][nu]=trace(Ap[m][mu]*Ap[phases.numMom()-m][nu]);
	}

	if(Layout::primaryNode())
	{
		general_data_base io_ap;
		sprintf(io_ap.name,"%s",(params.filename+"_ap.iog").c_str());
		QDPIO::cout << "write file to" << params.filename << std::endl;
		QDPIO::cout << "max q^2=" << params.max_mom2 << std::endl;
		io_ap.add_dimension(dim_momentum, mom_serial.size(),mom_serial.data());
		io_ap.add_dimension(dim_direction, Nd);
		io_ap.add_dimension(dim_temporary, 9);
		io_ap.add_dimension(dim_complex, 2);
		io_ap.initialize();
		int data_index=0;
		for(int m(0);m<mom_serial.size();m++)
		for(int dir=0; dir< Nd; ++dir)
		for(int ic2=0; ic2!=Nc; ic2++)
		for(int ic1=0; ic1!=Nc; ic1++)
		{
			io_ap.data[data_index]=Ap[m][dir].elem().elem().elem(ic1,ic2).real();
			data_index++;
			io_ap.data[data_index]=Ap[m][dir].elem().elem().elem(ic1,ic2).imag();
			data_index++;
		}
		QDPIO::cout << "writing A(p) to" << io_ap.name << std::endl;
		io_ap.save();
		
		general_data_base io_prop;
		sprintf(io_prop.name,"%s",(params.filename+"_prop.iog").c_str());
		QDPIO::cout << "write file to" << params.filename << std::endl;
		QDPIO::cout << "max q^2=" << params.max_mom2 << std::endl;
		io_prop.add_dimension(dim_momentum, (1+mom_serial.size())/2,mom_serial.data());
		io_prop.add_dimension(dim_direction, Nd);
		io_prop.add_dimension(dim_temporary, Nd);
		io_prop.add_dimension(dim_complex, 2);
		io_prop.initialize();
		data_index=0;
		for(int m(0); m<(1+mom_serial.size())/2;m++)
		for(int dir=0; dir< Nd; ++dir)
		for(int dir2=0; dir2<Nd; ++dir2)
		{
			io_prop.data[data_index]=Dp[m][dir][dir2].elem().elem().elem().real();
			data_index++;
			io_prop.data[data_index]=Dp[m][dir][dir2].elem().elem().elem().imag();
			data_index++;
		}
		QDPIO::cout << "writing file to" << io_prop.name << std::endl;
		io_prop.save();
	}
   
    
    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    
              
    pop(xml_out);  // npr

    snoop.stop();
    QDPIO::cout << InlineGluonNprPDFEnv::name << ": total time = "
		<< snoop.getTimeInSeconds()
		<< " secs" << std::endl;
    
    QDPIO::cout << InlineGluonNprPDFEnv::name << ": ran successfully" << std::endl;
    
    END_CODE();
  }
  
} //name space Chroma
