.PHONY:	run, pull, help

TARGETDIR 	:= ./pmopt
TARGET		:= ./pmopt
BRANCH		:= main

GNUMAKEFLAGS 	:= --no-print-directory

run:
	nohup sh run_all.sh &


pull:
	cd $(TARGETDIR) && git pull origin $(BRANCH) && make MODE=Release


help:
	$(TARGETDIR)/$(TARGET)