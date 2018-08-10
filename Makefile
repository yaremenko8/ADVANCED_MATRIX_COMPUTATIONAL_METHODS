CC = ghc
HSSDIRS =
HSSFLAGS =
all : mtxprog
mtxprog : main.c CMatrix.o Matrix.o FFI.o
	$(CC) --make -no-hs-main -optc-O $^ -o $@
CMatrix.o : CMatrix.hs Matrix.o FFI.o
	$(CC) $(HSSDIRS) $(HSSFLAGS) -c -O $<
Matrix.o : Matrix.hs
	$(CC) $(HSSDIRS) $(HSSFLAGS) -c -O $<
FFI.o : FFI.hs
	$(CC) $(HSSDIRS) $(HSSFLAGS) -c -O $<
clean :
	rm -f *.hi *.o CMatrix_stub.h mtxprog
