EXEBrute = knn-bruteforce
EXEDriver = knn-driver
runB: $(EXEBrute)
	futhark dataset -b -g 1i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk1d15.txt
	futhark dataset -b -g 60i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk60d15.txt
	futhark dataset -b -g 120i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk120d15.txt
	futhark dataset -b -g 1i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk1d50.txt
	futhark dataset -b -g 60i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk60d50.txt
	futhark dataset -b -g 120i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk120d50.txt
	futhark dataset -b -g 1i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk1d150.txt
	futhark dataset -b -g 60i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk60d150.txt
	futhark dataset -b -g 120i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> resk120d150.txt

runD: $(EXEDriver)
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5 < tmpdata  >> t1resk1d15.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5 < tmpdata >> t1resk60d15.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  -r 5 < tmpdata >> t1resk120d15.txt
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5  < tmpdata >> t1resk1d50.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5  < tmpdata >> t1resk60d50.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5  < tmpdata >> t1resk120d50.txt
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5  < tmpdata >> t1resk1d150.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5  < tmpdata >> t1resk60d150.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 5  < tmpdata >> t1resk120d150.txt

runDt10: $(EXEDriver)
	futhark dataset -b -g 10i32 -g 1i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n < tmpdata >> t10resk1d15.txt
	futhark dataset -b -g 10i32 -g 60i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n < tmpdata >> t10resk60d15.txt
	futhark dataset -b -g 10i32 -g 120i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n < tmpdata >> t10resk120d15.txt
	futhark dataset -b -g 10i32 -g 1i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  < tmpdata >> t10resk1d50.txt
	futhark dataset -b -g 10i32 -g 60i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  < tmpdata >> t10resk60d50.txt
	futhark dataset -b -g 10i32 -g 120i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  < tmpdata >> t10resk120d50.txt
	futhark dataset -b -g 10i32 -g 1i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  < tmpdata >> t10resk1d150.txt
	futhark dataset -b -g 10i32 -g 60i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  < tmpdata >> t10resk60d150.txt
	futhark dataset -b -g 10i32 -g 120i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  < tmpdata >> t10resk120d150.txt


#/$(EXEBrute) -t /dev/stdout -n -r 10 < tmpdata >> resk1d50.txt
# --runs=NUM
clean: 
	rm -f *.c
	rm -f *.actual
	rm -f *.expected
