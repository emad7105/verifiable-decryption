.PHONY: all clean

# all: vdec_params.h
all: vdec_params_tbox.h

# _fhe_example_params_cffi.o: fhe_example_params.h ../_lazer_cffi.o
#	python3 ../params_cffi_build.py fhe_example_params.h ../..

vdec_params_tbox.h: #fhe_example_params.py
	cd ../scripts && sage lnp-tbox-codegen.sage ../vdec/vdec_params_tbox.sage > ../vdec/vdec_params_tbox.h
#	cd ../scripts && sage lin-codegen.sage ../vdec/vdec_params.sage > ../vdec/vdec_params.h

#../_lazer_cffi.o:
#	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
