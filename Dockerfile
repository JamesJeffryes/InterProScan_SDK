FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update

# -----------------------------------------

RUN cpanm -i Config::IniFiles
RUN apt-get install libffi-dev libssl-dev
RUN pip install --upgrade requests[security]
RUN apt-add-repository ppa:webupd8team/java
RUN apt-get update
RUN apt-get -q install -y oracle-java8-installer
RUN apt-get install oracle-java8-set-default
RUN java -version
RUN which java

# Copy local wrapper files, and build

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make

RUN rm -rf /kb/runtime/java

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]