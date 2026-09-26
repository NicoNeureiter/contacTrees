# Dockerfile to build container for unit testing

FROM debian:stable

RUN apt-get update \
	&& apt-get install -y openjdk-25-jdk openjfx maven \
	&& rm -rf /var/lib/apt/lists/*

WORKDIR /root

COPY . ./contacTrees
WORKDIR /root/contacTrees

ENTRYPOINT ["mvn", "test"]