FROM broadinstitute/gtex_eqtl:V8

RUN cd /opt \
  && wget https://github.com/alyenkin/fastqtl/archive/master.zip \
  && unzip master.zip \
  && mv fastqtl-master fastqtl-sv \
  && cd fastqtl-sv \
  && mkdir obj \
  && mkdir bin \
  && make cleanall \
  && make RMATH=$RMATH
