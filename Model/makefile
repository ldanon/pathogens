 #Dengue.o: *.cpp *.h
 #	g++ -o Dengue.o  *.cpp *.h -I/opt/local/include/ -lm -L/opt/local/lib -lgsl -O5

a.out: *.cpp *.h
	g++ *.cpp  -lm -L/usr/local/lib -lgsl -lgslcblas

clean: 
	rm a.out
