ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

SOURCE_POSTFIX=c

INSTALL_DIR?=$(ROOT_DIR)/install
INCLUDE_DIR=$(ROOT_DIR)/include
SRC_DIR=$(ROOT_DIR)/src
OBJ_DIR=$(ROOT_DIR)/obj

EXE_NAME=fmax

CC=gcc
CFLAGS=-I$(INCLUDE_DIR) -L$(LIB_DIR) -std=c11 -O3 -Wall
CLIBS=-lm -lz
CINPUT=$(SRC_DIR)/*.$(SOURCE_POSTFIX)



all: compile link



compile:
	@cd $(OBJ_DIR) && \
	$(CC) $(CFLAGS) -c $(CINPUT) && \
	cd $(ROOT_DIR)

link:
	@$(CC) $(CFLAGS) $(OBJ_DIR)/* -o $(INSTALL_DIR)/$(EXE_NAME) $(CLIBS)

build:
	@$(CC) $(CFLAGS) $(CINPUT) -o $(INSTALL_DIR)/$(EXE_NAME) $(CLIBS)



clean_obj:
	@rm -f $(OBJ_DIR)/*

clean_install:
	@rm -rf $(INSTALL_DIR)/*

clean: clean_obj

clean_all: clean_obj clean_install
