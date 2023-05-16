EXEBrute = knn-bruteforce
EXEDriver = knn-driver
runB: $(EXEBrute)
	futhark dataset -b -g 1i64 -g [10][15]f32 -g [50][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk1d50.txt
	futhark dataset -b -g 3i64 -g [10][10]f32 -g [10][10]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk15d50.txt
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk60d50.txt
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk120d50.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk1d100.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk15d100.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk60d100.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk100d100.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk1d150.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk15d150.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk60d150.txt
#	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk100d150.txt

runD: $(EXEDriver)
	futhark dataset -b -g 1i32 -g 1i64 -g [10][15]f32 -g [50][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n -r 10 < tmpdata >> resk1d50.txt
# --runs=NUM
clean: 
	rm -f *.c
	rm -f *.actual
	rm -f *.expected