.PHONY: image open cluster push_image pm_membership
CR := "docker"
IMAGE_NAME := "open_clusters"
IMAGE_TAG := "1.0.0"
IMAGE := $(IMAGE_NAME):$(IMAGE_TAG)

image:
	$(CR) build --platform linux/amd64 -t $(IMAGE) . #If you are building an image on your M1 Mac to be used on Linux then use the --platform linux/amd64 to build a version for Intel chips

run: image
	docker run --platform linux/amd64 -v ~/open_clusters/data:/root/open_clusters/data -it $(IMAGE)
