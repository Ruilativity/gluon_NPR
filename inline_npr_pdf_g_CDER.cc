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
		
	if (paramtop.count("lmax") == 1) read(paramtop, "lmax", lmax);
	else lmax=-1;

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
    QDP::write(xml_out, "lmax", lmax);
    QDP::write(xml_out, "max_mom2", max_mom2);
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }

template <typemane T> multi4d<T> fast_FT(multi4d<T> to_FT)
{
	const Real twopi = 6.283185307179586476925286;
	multi4d<T> FT_result;
	FT_result.resize(to_FT.size4,to_FT.size3,to_FT.size2,to_FT.size1)
	multi4d<T> to_FT_sub1, to_FT_sub2, to_FT_sub3;
	multi4d<T> FT_result_sub1, FT_result_sub2, FT_result_sub3;
	if (to_FT.size4 > 1){
		if (to_FT.size4%3 == 0){
			to_FT_sub1.resize(to_FT.size4/3,to_FT.size3,to_FT.size2,to_FT.size1)
			to_FT_sub2.resize(to_FT.size4/3,to_FT.size3,to_FT.size2,to_FT.size1)
			to_FT_sub3.resize(to_FT.size4/3,to_FT.size3,to_FT.size2,to_FT.size1)
			for(int i=0;i<to_FT_sub1.size4;i++)
			for(int j=0;j<to_FT_sub1.size3;j++)
			for(int k=0;k<to_FT_sub1.size2;k++)
			for(int l=0;l<to_FT_sub1.size1;l++){
				to_FT_sub1(i,j,k,l)=to_FT(3*i,j,k,l);
				to_FT_sub2(i,j,k,l)=to_FT(3*i+1,j,k,l);
				to_FT_sub3(i,j,k,l)=to_FT(3*i+2,j,k,l);
			}
			FT_result_sub1=fast_FT(to_FT_sub1);
			FT_result_sub2=fast_FT(to_FT_sub2);
			FT_result_sub3=fast_FT(to_FT_sub3);
			for(int i=0;i<FT_result_sub1.size4;i++)
			for(int j=0;j<FT_result_sub1.size3;j++)
			for(int k=0;k<FT_result_sub1.size2;k++)
			for(int l=0;l<FT_result_sub1.size1;l++){
				FT_result(i,j,k,l)=FT_result_sub1(i,j,k,l) + FT_result_sub2(i,j,k,l)*cmplx(cos(twopi*i/FT_result.size4),sin(twopi*i/FT_result.size4)) + FT_result_sub3(i,j,k,l)*cmplx(cos(2*twopi*i/FT_result.size4),sin(2*twopi*i/FT_result.size4));
				FT_result(i+FT_result_sub1.size4,j,k,l)=FT_result_sub1(i,j,k,l) + FT_result_sub2(i,j,k,l)*cmplx(cos(twopi*i/FT_result.size4+twopi/3),sin(twopi*i/FT_result.size4+twopi/3)) + FT_result_sub3(i,j,k,l)*cmplx(cos(2*twopi*i/FT_result.size4+2*twopi/3),sin(2*twopi*i/FT_result.size4+2*twopi/3));
				FT_result(i+2*FT_result_sub1.size4,j,k,l)=FT_result_sub1(i,j,k,l) + FT_result_sub2(i,j,k,l)*cmplx(cos(twopi*i/FT_result.size4+2*twopi/3),sin(twopi*i/FT_result.size4+2*twopi/3)) + FT_result_sub3(i,j,k,l)*cmplx(cos(2*twopi*i/FT_result.size4+twopi/3),sin(2*twopi*i/FT_result.size4+twopi/3));
			}
			
		}
		else if(to_FT.size4%2 == 0){
			to_FT_sub1.resize(to_FT.size4/2,to_FT.size3,to_FT.size2,to_FT.size1)
			to_FT_sub2.resize(to_FT.size4/2,to_FT.size3,to_FT.size2,to_FT.size1)
			for(int i=0;i<to_FT_sub1.size4;i++)
			for(int j=0;j<to_FT_sub1.size3;j++)
			for(int k=0;k<to_FT_sub1.size2;k++)
			for(int l=0;l<to_FT_sub1.size1;l++){
				to_FT_sub1(i,j,k,l)=to_FT(2*i,j,k,l);
				to_FT_sub2(i,j,k,l)=to_FT(2*i+1,j,k,l);
			}
			FT_result_sub1=fast_FT(to_FT_sub1);
			FT_result_sub2=fast_FT(to_FT_sub2);
			for(int i=0;i<FT_result_sub1.size4;i++)
			for(int j=0;j<FT_result_sub1.size3;j++)
			for(int k=0;k<FT_result_sub1.size2;k++)
			for(int l=0;l<FT_result_sub1.size1;l++){
				FT_result(i,j,k,l)=FT_result_sub1(i,j,k,l) + FT_result_sub2(i,j,k,l)*cmplx(cos(twopi*i/FT_result.size4),sin(twopi*i/FT_result.size4)) ;
				FT_result(i+FT_result_sub1.size4,j,k,l)=FT_result_sub1(i,j,k,l) - FT_result_sub2(i,j,k,l)*cmplx(cos(twopi*i/FT_result.size4),sin(twopi*i/FT_result.size4));
			}
			
		}
		
	}
}

template <typemane T> multi4d<T> fast_FT(multi4d<T> to_FT, int sign_flag) // sign_flag=1 for FT, -1 for inverse FT
{
	const Real twopi = 6.283185307179586476925286;
	multi4d<T> FT_result;
	FT_result.resize(to_FT.size4,to_FT.size3,to_FT.size2,to_FT.size1);
	FT_result=0;
	multi1d<multi4d<T>> to_FT_sub;
	multi1d<multi4d<T>> FT_result_sub;
	multi1d<int> sizes(4);
	multi1d<int> sub_indices(4);
	multi1d<int> sub_in_full_indices(4);
	sizes[0]=to_FT.size4;
	sizes[1]=to_FT.size3;
	sizes[2]=to_FT.size2;
	sizes[3]=to_FT.size1;
	for (int idir=0;idir<4;idir++){
		if (sizes[idir] > 1){
			for (int batches=2;batches<=3;batches++)
				if (sizes[idir]%batches == 0){
					to_FT_sub.resize[batches];
					FT_result_sub.resize[batches];
					multi1d<int> batch_sizes(sizes);
					batch_sizes(idir)/=batches;
					for (int ibatch=0;ibatch<batches;ibatch++){
						to_FT_sub[ibatch].resize(batch_sizes[0],batch_sizes[1],batch_sizes[2],batch_sizes[3]);
						for(int i=0;i<batch_sizes[0];i++){
							sub_indices[0]=i;sub_in_full_indices[0]=i;
							for(int j=0;j<batch_sizes[1];j++){
								sub_indices[1]=j;sub_in_full_indices[1]=j;
								for(int k=0;k<batch_sizes[2];k++){
									sub_indices[2]=k;sub_in_full_indices[2]=k;
									for(int l=0;l<batch_sizes[3];l++){
										sub_indices[3]=l;sub_in_full_indices[3]=l;
										sub_in_full_indices[idir]=batches*sub_indices[idir]+ibatch;
										to_FT_sub[ibatch](sub_indices[0],sub_indices[1],sub_indices[2],sub_indices[3])=to_FT(sub_in_full_indices[0],sub_in_full_indices[1],sub_in_full_indices[2],sub_in_full_indices[3]);
									}
								}
							}
						}
						FT_result_sub[ibatch]=fast_FT(to_FT_sub[ibatch]);
					}
					for (int ibatch1=0;ibatch1<batches;ibatch1++){
						for(int i=0;i<batch_sizes[0];i++){
							sub_indices[0]=i;sub_in_full_indices[0]=i;
							for(int j=0;j<batch_sizes[1];j++){
								sub_indices[1]=j;sub_in_full_indices[1]=j;
								for(int k=0;k<batch_sizes[2];k++){
									sub_indices[2]=k;sub_in_full_indices[2]=k;
									for(int l=0;l<batch_sizes[3];l++){
										sub_indices[3]=l;sub_in_full_indices[3]=l;
										sub_in_full_indices[idir]+=ibatch1*batch_sizes[idir];
										for (int ibatch2=0;ibatch2<batches;ibatch2++){
											FT_result(sub_in_full_indices[0],sub_in_full_indices[1],sub_in_full_indices[2],sub_in_full_indices[3])+=FT_result_sub[ibatch2](sub_indices[0],sub_indices[1],sub_indices[2],sub_indices[3]) *cmplx(cos(twopi*sub_indices[ibatch2]/sizes[ibatch2]+ibatch*twopi/batches),sign_flag*sin(twopi*sub_indices[ibatch2]/sizes[ibatch2]+ibatch*twopi/batches));
										}
									}
								}
							}
						}
					}
					break;
				}
			
		}
		else {
			if(idir<3) continue;
			else FT_result=to_FT;
		}
	}
	return FT_result;
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
    SftMom phases(params.max_mom2, false,  Nd); // 4D fourier transform

	std::vector<int> mom_serial(phases.numMom());

	  double g0=1.0;
	  LatticeColorMatrix trUmat;
	  multi1d<LatticeColorMatrix> ai(Nd);
	  multi1d<multi4d<PColorMatrix< RComplex<REAL>, Nc>>> ai_x(Nd);
	  for(int mu=0; mu<Nd; ++mu){
		  trUmat=2.0/Nc*trace(u[mu]);
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

//for(int mu=0; mu<Nd; ++mu) checkpoint[mu]=1;
	  
	  

	multi2d<ColorMatrix> Ap(phases.numMom(),Nd);
	multi4d<double> Dp((phases.numMom()+1)/2,Nd ,Nd,2);
	Complex shift_phase;
	for (int m=0; m < phases.numMom(); m++){
		mom_serial[m]=(50-phases.numToMom(m)[0])+(50-phases.numToMom(m)[1])*100+(50-phases.numToMom(m)[2]) *10000+(50-phases.numToMom(m)[3])*1000000;
		
		for(int mu=0; mu<Nd; ++mu){
			const Real twopi = 6.283185307179586476925286;
			Real p_dot_x;
			p_dot_x=phases.numToMom(m)[mu]*twopi/Layout::lattSize()[mu]/2.0;
			shift_phase=cmplx(cos(p_dot_x),sin(p_dot_x));
			Ap[m][mu] = shift_phase*sum(phases[m]*ai[mu]);
		}
	}
	for (int m=0; m < (1+phases.numMom())/2; m++){
		for(int mu=0; mu<Nd; ++mu)
			for(int nu=0; nu<Nd; ++nu){
				Dp[m][mu][nu][0]=0.;
				Dp[m][mu][nu][1]=0.;
				for(int i=0; i<Nc; ++i)
					for(int j=0; j<Nc; ++j){
						Dp[m][mu][nu][0]+=(Ap[m][mu].elem().elem().elem(i,j).real()*Ap[phases.numMom()-m-1][nu].elem().elem().elem(j,i).real()-Ap[m][mu].elem().elem().elem(i,j).imag()*Ap[phases.numMom()-m-1][nu].elem().elem().elem(j,i).imag());
                                		Dp[m][mu][nu][1]+=(Ap[m][mu].elem().elem().elem(i,j).imag()*Ap[phases.numMom()-m-1][nu].elem().elem().elem(j,i).real()+Ap[m][mu].elem().elem().elem(i,j).real()*Ap[phases.numMom()-m-1][nu].elem().elem().elem(j,i).imag());
					}
			}
	}

	  
	  Double w_plaq, s_plaq, t_plaq, link;
      multi2d<Double> plane_plaq;
	  MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
	  multi2d<LatticeColorMatrix> plane_plaq_matrix;
	  multi1d<multi4d<Double>> gluon_o;
	  multi1d<LatticeReal> gluon_o2,gluon_o3,gluon_o4;
	  
	  if(params.lmax>=0){
	  plane_plaq_matrix.resize(Nd,Nd);
	  
	  for(int mu=1; mu < Nd; ++mu)
	  {
		  for(int nu=0; nu < mu; ++nu)
		  {
			  plane_plaq_matrix[mu][nu]=u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]);
			  plane_plaq_matrix[nu][mu]=adj(plane_plaq_matrix[mu][nu]);
		  }
	  }
	  multi2d<LatticeColorMatrix> FF,FF00,FF01,FF10,FF11;
	  FF.resize(Nd,Nd);
	  FF00.resize(Nd,Nd);
	  FF01.resize(Nd,Nd);
	  FF10.resize(Nd,Nd);
	  FF11.resize(Nd,Nd);
	  for(int mu=1; mu < Nd; ++mu)
	  {
		  FF[mu][mu]=0;
		  for(int nu=0; nu < mu; ++nu)
		  {
			  FF00[mu][nu]=(plane_plaq_matrix[mu][nu]-plane_plaq_matrix[nu][mu])/2.0;
			  FF01[mu][nu]=shift(FF00[mu][nu],BACKWARD,mu);
			  FF10[mu][nu]=shift(FF00[mu][nu],BACKWARD,nu);
			  FF11[mu][nu]=shift(FF01[mu][nu],BACKWARD,nu);
			  FF[mu][nu]=(FF00[mu][nu]+FF01[mu][nu]+FF10[mu][nu]+FF11[mu][nu])/4.0;
			  FF[nu][nu]=0;
		  }
	  }
	  
	  gluon_o1.resize(params.lmax);
	  gluon_o2.resize(params.lmax);
	  gluon_o3.resize(params.lmax);
	  gluon_o4.resize(params.lmax);
	  multi2d<LatticeColorMatrix> FF_shift;
	  FF_shift.resize(Nd,Nd);
	  LatticeColorMatrix tmp;
	  LatticeColorMatrix u_shift;
	  u_shift=1;
	  for(int mu=0; mu<Nd; ++mu)
		   for(int nu=0; nu<mu; ++nu)
		   {
			  FF_shift[mu][nu]=FF[mu][nu];
		   }
	  for(int l=0; l < params.lmax; ++l)
	  {
		  gluon_o[l]=0;gluon_o1[l]=0;gluon_o2[l]=0;gluon_o3[l]=0;gluon_o4[l]=0;
		  for(int mu=0; mu<Nd; ++mu)
		  {
			  if(mu != 3)
			  {
				  gluon_o3[l]+=real(trace(FF[mu][3]*u_shift*FF_shift[mu][3]*adj(u_shift)));
			  }
			  for(int nu=0; nu<Nd; ++nu)
				 if(nu != mu) gluon_o4[l]+=real(trace(FF[mu][nu]*u_shift*FF_shift[mu][nu]*adj(u_shift)))/Nd;
		  }
		  gluon_o2[l]=gluon_o3[l]-gluon_o4[l];
		  gluon_o[l].resize(Layout::lattSize()[0],Layout::lattSize()[1],Layout::lattSize()[2],Layout::lattSize()[3]);
		  for(int i=0; i<Layout::lattSize()[0];++i)
		  for(int j=0; j<Layout::lattSize()[1];++j)
		  for(int k=0; k<Layout::lattSize()[2];++k)
		  for(int t=0; t<Layout::lattSize()[3];++t){
			  coord[0]=i;
			  coord[1]=j;
			  coord[2]=k;
			  coord[3]=t;
			  gluon_o[l](i,j,k,t)=peekSite(gluon_o2[l],coord);
		  }
		  gluon_o[l]=gluon_o1[l]-gluon_o2[l];
		  for(int mu=0;mu<Nd;++mu)
			  for(int nu=0; nu<mu; ++nu){
				  tmp=shift(FF_shift[mu][nu],FORWARD,2);
				  FF_shift[mu][nu]=tmp;
			  }
		  if(l>0) tmp=u[2]*shift(u_shift,FORWARD,2);
		  else tmp=u[2];
		  u_shift=tmp;
		  
	  }
		  std::string out_fname_c("data_test_op.dat");
		  TextFileWriter fout(out_fname_c);
		  fout << "#t z O_tt O_mumu/4 O_tt-O_mumu/4" << "\n";
		  for(int l=0; l < params.lmax; ++l)
		  for(int t=0; t<Layout::lattSize()[3];++t)
		  {
				  fout <<t<<"\t"<<l<<"\t"<< gluon_o3_t[t][l]<<"\t"<< gluon_o4_t[t][l]<<"\t"<< gluon_o3_t[t][l]-gluon_o4_t[t][l]<< "\n";
		  }
		  fout.close();
	  }
	  
	  
	if(Layout::primaryNode())
	{
		general_data_base io_ap;
		sprintf(io_ap.name,"%s",(params.filename+"_ap.iog").c_str());
		QDPIO::cout << "write file to" << io_ap.name << std::endl;
		QDPIO::cout << "max q^2=" << params.max_mom2 << std::endl;
		io_ap.add_dimension(dim_momentum, mom_serial.size(),mom_serial.data());
		io_ap.add_dimension(dim_direction, Nd);
		io_ap.add_dimension(dim_temporary, 9);
		io_ap.add_dimension(dim_complex, 2);
		io_ap.initialize();
		int data_index=0;
		for(int m(0);m<mom_serial.size();m++)
		for(int dir=0; dir< Nd; ++dir)
		for(int ic1=0; ic1!=Nc; ic1++)
		for(int ic2=0; ic2!=Nc; ic2++)
		{
			io_ap.data[data_index]=Ap[m][dir].elem().elem().elem(ic1,ic2).real();
			data_index++;
			io_ap.data[data_index]=Ap[m][dir].elem().elem().elem(ic1,ic2).imag();
			data_index++;
		}
		QDPIO::cout << "writing A(p) to" << io_ap.name << std::endl;
		io_ap.save();
		
		if(params.lmax>=0){
		general_data_base io_op;
		sprintf(io_op.name,"%s",(params.filename+"_op.iog").c_str());
		QDPIO::cout << "write file to" << io_op.name << std::endl;
		QDPIO::cout << "max link=" << params.lmax << std::endl;
		io_op.add_dimension(dim_operator, 3);
		io_op.add_dimension(dim_displacement, params.lmax);
		io_op.initialize();
		data_index=0;
		for(int l(0); l<params.lmax; l++)
		{
			io_op.data[data_index]=gluon_o[l].elem().elem().elem().elem();
			data_index++;
		}
		for(int l(0); l<params.lmax; l++)
		{
			io_op.data[data_index]=gluon_o1[l].elem().elem().elem().elem();
			data_index++;
		}
		for(int l(0); l<params.lmax; l++)
		{
			io_op.data[data_index]=gluon_o2[l].elem().elem().elem().elem();
			data_index++;
		}
		QDPIO::cout << "writing file to" << io_op.name << std::endl;
		io_op.save();
		}
		
		general_data_base io_prop;
		sprintf(io_prop.name,"%s",(params.filename+"_prop.iog").c_str());
		QDPIO::cout << "write file to" << io_prop.name << std::endl;
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
			io_prop.data[data_index]=Dp[m][dir][dir2][0];
			data_index++;
			io_prop.data[data_index]=Dp[m][dir][dir2][1];
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
