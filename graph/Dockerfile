# base this container on neo4j version 3.5.4
FROM neo4j:3.5.4

# get some credit
LABEL maintainer="powen@renci.org"

# include java version 11
RUN apk --no-cache add openjdk11 --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community

# Set the install path to APOC for neo4j
ENV APOC_URI https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/3.5.0.2/apoc-3.5.0.2-all.jar

# create a directory for the ssl certificate
RUN mkdir -p /ssl

# make sure it is accessable
RUN chown neo4j /ssl

# make a directory for the APOC plugin
RUN mkdir /plugins

# get the APOC plugin for neo4j
RUN wget $APOC_URI \
    && mv apoc-3.5.0.2-all.jar /plugins

# expose the common neo4j ports
EXPOSE 7474 7473 7687

# start up neo4j on the container
CMD ["neo4j"]