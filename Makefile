TOPDIR=${WORK}/boram/bin2
#CHROMA=/ccs/home/ybyang1/scratch/chroma_single_devel
#CHROMA=/global/project/projectdirs/m2838/chroma_single
#CHROMA=${TOPDIR}/chroma_LANL-double
#CHROMA=${TOPDIR}/pdfchroma
CHROMA=${TOPDIR}/chroma-double
CONFIG=$(CHROMA)/bin/chroma-config

MGLIBS= -lwilsonmg -lqopqdp -lqdp_common -lqdp_int -lqdp_f -lqdp_d -lqdp_df -lqdp_f2 -lqdp_d2 -lqdp_df2 -lqdp_f3 -lqdp_d3 -lqdp_df3 -lqdp_fn -lqdp_dn -lqdp_dfn -lqio -llime -lqla_c99 -lqla_cmath -lqla_d2 -lqla_d3 -lqla_d -lqla_df2 -lqla_df3 -lqla_df -lqla_dfn -lqla_dn -lqla_dq2 -lqla_dq3 -lqla_dq -lqla_dqn -lqla_f2 -lqla_f3 -lqla_f -lqla_fn -lqla_int -lqla_q2 -lqla_q3 -lqla_q -lqla_qn -lqla_random
MGLDFLAGS=-L$(TOPDIR)/qla/lib -L$(TOPDIR)/lib -L$(TOPDIR)/qopqdp/lib
MGCXXFLAGS=-L$(TOPDIR)/qla/lib -L$(TOPDIR)/qdpc/lib -L$(TOPDIR)/qopqdp/lib

CXX=$(shell $(CONFIG) --cxx) 
#CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I$(TOPDIR)/qla/include -I. $(MGCXXFLAGS) 
CXXFLAGS=-D_PDF_ $(shell $(CONFIG) --cxxflags) -I$(TOPDIR)/qla/include -I. $(MGCXXFLAGS) 
LDFLAGS=$(shell $(CONFIG) --ldflags) $(MGLDFLAGS)
LIBS=$(shell $(CONFIG) --libs) $(MGLIBS)

HDRS=inline_npr_pdf_g.h \


OBJS=chroma.o \
     io_general.o io_general_class.o\
     inline_npr_pdf_g.o \

chroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)


%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -rf chroma $(OBJS) *~
