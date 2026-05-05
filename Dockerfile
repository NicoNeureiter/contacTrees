# Dockerfile to build container for unit testing

FROM debian:stable

RUN apt-get update \
	&& apt-get install -y git openjdk-25-jdk openjfx maven \
	&& rm -rf /var/lib/apt/lists/*

WORKDIR /root

ARG FEAST_REF="beast2.8-migration"
RUN git clone https://github.com/tgvaughan/feast.git /root/feast \
    && cd /root/feast \
    && git checkout "$FEAST_REF" \
    && mvn install -DskipTests

COPY . ./contacTrees
WORKDIR /root/contacTrees

ENTRYPOINT ["mvn", "test"]