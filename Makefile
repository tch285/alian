REPO  ?= tch285
IMAGE ?= alian
TAG   ?= latest

build:
	docker build -t $(REPO)/$(IMAGE):$(TAG) .
push:
	docker push $(REPO)/$(IMAGE):$(TAG)
run:
	docker run -it --rm $(REPO)/$(IMAGE):$(TAG) /bin/bash -l

.PHONY: build push run