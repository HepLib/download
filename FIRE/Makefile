ifeq ($(strip $(MALLOC)),)
MALLOC:=-ljemalloc
endif

# Note: one needs to set INSTALL_PATH
export MALLOC

defalut: oM oF oQ oP

oM:
	$(MAKE) -C M
oF:
	$(MAKE) -C F
oQ:
	$(MAKE) -C Q
oP:
	$(MAKE) -C P
		
cleanall: clean cleanexe

clean:
	$(MAKE) clean -C M
	$(MAKE) clean -C F
	$(MAKE) clean -C Q
	$(MAKE) clean -C P

cleanexe:
	rm -f M/FIRE
	rm -f F/FIRE
	rm -f Q/FIRE
	rm -f P/FIRE

