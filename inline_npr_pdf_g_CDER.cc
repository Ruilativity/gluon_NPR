/*! \file Id: inline_npr_pdf_w.cc  Modified by R. Zhang in Nov. 1st, 2019
 * \brief Inline construction of NPR propagator
 *
 * Propagator calculations
 */
#include <fftw3.h>
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

//convert the lattice color matrix to fftw arrays, waiting for FT.
multi1d<fftw_complex*> lattice_color_matrix_to_fftw_array(LatticeColorMatrix ai_x){
	multi1d<multi1d<Complex>> color_array(Nc*Nc);
	multi1d<fftw_complex*> color_fftw_array(Nc*Nc);
	int color_index=0;
	for(int i=0;i<Nc;i++){
		for(int j=0;j<Nc;j++){
			color_array[color_index].resize(Layout::vol());
			QDP_extract(color_array[color_index],peekColor(ai_x,i,j),all); // order: x(fastest),y,z,t(slowest)
			color_fftw_array[color_index] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
			for(int k=0;k<color_array[color_index].size();k++){
				color_fftw_array[color_index][k][0]=color_array[color_index][k].real();
				color_fftw_array[color_index][k][1]=color_array[color_index][k].imag();
			}
			color_index++;
		}
	}
	return color_fftw_array;
}
//convert the array of lattice complex to fftw arrays, waiting for FT.
multi1d<fftw_complex*> lattice_complex_to_fftw_array(multi1d<LatticeComplex> op){
	multi1d<multi1d<Complex>> op_array;
	multi1d<fftw_complex*> op_fftw_array;
	for(int i=0;i<op_array.size();i++){
		op_array[i].resize(Layout::vol());
		QDP_extract(op_array[i],op[i],all); // order: x(fastest),y,z,t(slowest)
		op_fftw_array[i] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
		for(int k=0;k<op_array[i].size();k++){
			op_fftw_array[i][k][0]=op_array[i][k].real();
			op_fftw_array[i][k][1]=op_array[i][k].imag();
		}
	}
	return op_fftw_array;
}
//convert the lattice complex to fftw arrays, waiting for FT.
fftw_complex* lattice_complex_to_fftw_array(LatticeComplex op){
	multi1d<Complex> op_array;
	fftw_complex* op_fftw_array;
	op_array.resize(Layout::vol());
	QDP_extract(op_array,op,all); // order: x(fastest),y,z,t(slowest)
	op_fftw_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
	for(int k=0;k<op_array.size();k++){
		op_fftw_array[k][0]=op_array[k].real();
		op_fftw_array[k][1]=op_array[k].imag();
	}
	return op_fftw_array;
}

// FT (sign=1) or inverse FT (sign=-1), replace the input array with its FT. (normalize the data with 1/V after inverse FT)
void fftw_transformation(multi1d<fftw_complex*> color_fftw_in_array, int sign, bool norm_flag){
	multi1d<fftw_complex*> color_fftw_out_array(color_fftw_in_array.size());
	for(int i=0;i<color_fftw_in_array.size();i++){
		color_fftw_out_array[i] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
	}
	int sizes[Nd];
	for(int i=0;i<Nd;i++) sizes[i]=Layout::lattSize()[i]; // order: x(fastest),y,z,t(slowest)
	for(int i=0;i<color_fftw_in_array.size();i++){
		fftw_plan p=fftw_plan_dft(Nd, sizes, color_fftw_in_array[i], color_fftw_out_array[i],
								sign, FFTW_MEASURE);
		fftw_execute(p);
		if(norm_flag) for(int j=0;j<sizeof(color_fftw_in_array[i]);j++) color_fftw_in_array[i]/=Layout::vol();
	}
	return color_fftw_out_array;
}
// FT (sign=1) or inverse FT (sign=-1), replace the input array with its FT. (normalize the data with 1/V after inverse FT)
void fftw_transformation(fftw_complex* scalar_fftw_in_array, int sign, bool norm_flag){
	fftw_complex* scalar_fftw_out_array;
	color_fftw_out_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
	int sizes[Nd];
	for(int i=0;i<Nd;i++) sizes[i]=Layout::lattSize()[i]; // order: x(fastest),y,z,t(slowest)
	fftw_plan p=fftw_plan_dft(Nd, sizes, scalar_fftw_in_array, scalar_fftw_out_array,
							sign, FFTW_MEASURE);
	fftw_execute(p);
	if(norm_flag) for(int j=0;j<sizeof(color_fftw_in_array);j++) color_fftw_in_array/=Layout::vol();
	return scalar_fftw_out_array;
}


/* CDER for LaMET operator.
	input: l links of the operator (lattice Complex), gauge field (4 direction lattice color matrix), p*4 momentum sets, r1 (A to O), r2 (A to A).
	output: l links * p momentums of the operator (Complex)
 */
multi1d<multi1d<Complex>> CDER_combination(multi1d<LatticeComplex> gluon_op_x, multi1d<LatticeColorMatrix> ai_x, multi1d<multi1d<int>> mom_sets, int r1, int r2){
	LatticeReal filter_r1, filter_r2;
	LatticeBoolean mask = false;
	//create filter |x|<r1: f(x)=\theta(r1-|x|)
	mask = (pow(Layout::latticeCoordinate(0),2)+pow(Layout::latticeCoordinate(1),2)+pow(Layout::latticeCoordinate(2),2)+pow(Layout::latticeCoordinate(3),2)) <= pow(r1,2);
	filter_r1=where(mask, 1.0, 0.0);

	//create filter |x|<r2: g(x)=\theta(r2-|x|)
	mask = (pow(Layout::latticeCoordinate(0),2)+pow(Layout::latticeCoordinate(1),2)+pow(Layout::latticeCoordinate(2),2)+pow(Layout::latticeCoordinate(3),2)) <= pow(r2,2);
	filter_r2=where(mask, 1.0, 0.0);
	
	//FT filter_r1 f(p)
	fftw_complex* filter_r1_fftw_x, filter_r1_fftw_p;
	filter_r1_fftw_x=lattice_complex_to_fftw_array(filter_r1);
	filter_r1_fftw_p=fftw_transformation(filter_r1_fftw_x, 1, false);
	
	//FT lamet operator O(p)
	multi1d<fftw_complex*> gluon_op_fftw_x, gluon_op_fftw_p;
	gluon_op_fftw_x = lattice_complex_to_fftw_array(gluon_op_x);
	gluon_op_fftw_p = fftw_transformation(gluon_op_fftw_x, 1, false);
	
	//multiply the two arrays, then inverse ft to coordinate space. \tilde{O}(x)=F[O(p)f(p)]
	multi1d<fftw_complex*> gluon_op_filtered_fftw_x(gluon_op_fftw_x.size()), gluon_op_filtered_fftw_p(gluon_op_fftw_x.size());
	
	for(int i=0; i < gluon_op_fftw_x.size();i++){
		gluon_op_filtered_fftw_x[i]= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
		for(int j=0; j< sizeof(filter_r1_fftw_x); j++){
			gluon_op_filtered_fftw_p[i][j][0]=gluon_op_fftw_p[i][j][0]*filter_r1_fftw_p[j][0]-gluon_op_fftw_p[i][j][1]*filter_r1_fftw_p[j][1];
			gluon_op_filtered_fftw_p[i][j][1]=gluon_op_fftw_p[i][j][0]*filter_r1_fftw_p[j][1]+gluon_op_fftw_p[i][j][1]*filter_r1_fftw_p[j][0];
		}
	}
	gluon_op_filtered_fftw_x=fftw_transformation(gluon_op_filtered_fftw_p, -1, true);
	
	
	//Convert gauge field A_\mu(x) to fftw array
	multi1d<multi1d<fftw_complex*>> ai_fftw_x(Nd), ai_fftw_p(Nd);
	for(int mu=0;mu<Nd;i++){
		ai_fftw_x[mu]=lattice_complex_to_fftw_array(ai_x[mu]);
		ai_fftw_p[mu]=fftw_transformation(ai_x[mu], 1, false);
	}
	
	/*
	 multiply the operator with gauge filed, then inverse ft to coordinate space. B(x)=\tilde{O}(x)A(x)
	 Dimensions: link length, spacetime direction of gauge filed, color index
	 */
	multi1d<multi1d<multi1d<fftw_complex*>>> gluon_op_filtered_ai_fftw_x(gluon_op_filtered_fftw_x.size()) gluon_op_filtered_ai_fftw_p(gluon_op_filtered_fftw_x.size());
	for(int i=0;i<gluon_op_filtered_fftw_x.size();i++){
		gluon_op_filtered_ai_fftw_x[i].resize(Nd);
		gluon_op_filtered_ai_fftw_p[i].resize(Nd);
		for(int j=0;j<Nd;j++){
			gluon_op_filtered_ai_fftw_x[i][j].resize(Nc*Nc);
			for(int k=0;k<Nc*Nc;k++){
				gluon_op_filtered_ai_fftw_x[i][j][k]=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Layout::vol());
				for(int site=0;site<Layout::vol();site++){
					gluon_op_filtered_ai_fftw_x[i][j][k][size][0]=gluon_op_filtered_fftw_x[i][size][0]*ai_fftw_x[j][k][size][0]-gluon_op_filtered_fftw_x[i][size][1]*ai_fftw_x[j][k][size][1];
					gluon_op_filtered_ai_fftw_x[i][j][k][size][1]=gluon_op_filtered_fftw_x[i][size][0]*ai_fftw_x[j][k][size][1]+gluon_op_filtered_fftw_x[i][size][1]*ai_fftw_x[j][k][size][0];
				}
			}
			gluon_op_filtered_ai_fftw_p[i][j]=fftw_transformation(gluon_op_filtered_ai_fftw_x[i][j],1,false);
		}
	}
	
	
	
	
	

}


//factorize one dimension and rearrange the coordinates for FFT.
std::vector<int> factorize(int length, std::vector<int> &index_list){
	int i=2;
	for(;i<=sqrt(length);i++){
		if(length%i ==0) {
			index_list.push_back(i);
			return factorize(length/i,index_list);
		}
	}
	if(i<= length)index_list.push_back(length);
	return index_list;
}
			
multi1d<int> re_arange(int length){
	multi1d<int> new_order;
	multi1d<int> tmp_order;
	new_order.resize(1);
	new_order=0;
	std::vector<int> factorization;
	factorization=factorize(length,factorization);
	for(int i=0;i<factorization.size();i++){
		length/=factorization[i];
		tmp_order=new_order;
		for(int j=1;j<factorization[i];j++)
			new_order=concat(new_order, tmp_order + (j-1)*length);
	}
	return new_order;
}

//fold 1-d array to multi-dimensional array

template <typemane T> multi1d<multi1d<T>> fold_1d(multi1d<T> to_fold, int degeneration){
	multi1d<multi1d<T>> folded;
	folded.resize(degeneration);
	int sub_len;
	if(to_fold.size()%degeneration != 0){
		QDPIO::cerr << InlineGluonNprPDFEnv::name << ": folding arrays of incorrect shape."
				<< std::endl;
		QDP_abort(1);
	}
	sublen=to_fold.size()/degeneration;
	int index=0;
	for(int i=0;i<degeneration;i++){
		folded[i].resize(sublen);
		for (int j=0; j<sublen;j++){
			folded[i][j]=to_fold[index];
			index++;
		}
	}
	return folded;
}


//FFT with arrangement for arrays, The data structure should be a multi1d array, can be nested.

/* not completed.
template <typemane T> multi1d<T> fast_FT_rearrange(multi1d<T> to_FT, int sign_flag) // sign_flag=1 for FT, -1 for inverse FT
{
	const Real twopi = 6.283185307179586476925286;
	multi1d<int> new_order=rearange(to_FT.size());
	
	return FT_result;
}
*/


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
