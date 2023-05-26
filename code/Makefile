EXEBrute = knn-bruteforce
EXEDriver = knn-driver
EXEAcc = knn-comparer
runB: $(EXEBrute)
	futhark dataset -b -g 5i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk5d15.txt
	futhark dataset -b -g 15i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk15d15.txt
	futhark dataset -b -g 30i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk30d15.txt
	futhark dataset -b -g 5i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk5d50.txt
	futhark dataset -b -g 15i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk15d50.txt
	futhark dataset -b -g 30i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk30d50.txt
	futhark dataset -b -g 5i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk5d150.txt
	futhark dataset -b -g 15i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk15d150.txt
	futhark dataset -b -g 30i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEBrute) -t /dev/stdout -n < tmpdata >> B_resk30d150.txt

test: $(EXEAcc) $(EXEDriver)
	futhark dataset -b -g 1i32 -g 5i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3 < tmpdata  >> t1resk5d15.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk5d15.txt

test2:
	echo 'hej'


accD: $(EXEAcc)
	#futhark dataset -b -g 1i32 -g 1i64 -g [500000][15]f32 -g [100][15]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk1d15.txt
	#futhark dataset -b -g 1i32 -g 60i64 -g [500000][15]f32 -g [100][15]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk60d15.txt
	#futhark dataset -b -g 1i32 -g 120i64 -g [500000][15]f32 -g [100][15]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk120d15.txt
	#futhark dataset -b -g 1i32 -g 5i64 -g [500000][15]f32 -g [100][15]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk5d15.txt
	#futhark dataset -b -g 1i32 -g 15i64 -g [500000][15]f32 -g [100][15]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk15d15.txt
	#futhark dataset -b -g 1i32 -g 30i64 -g [500000][15]f32 -g [100][15]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk30d15.txt
	#futhark dataset -b -g 1i32 -g 1i64 -g [500000][50]f32 -g [100][50]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk1d50.txt
	#futhark dataset -b -g 1i32 -g 5i64 -g [500000][50]f32 -g [100][50]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk5d50.txt
	#futhark dataset -b -g 1i32 -g 15i64 -g [500000][50]f32 -g [100][50]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk15d50.txt
	#futhark dataset -b -g 1i32 -g 30i64 -g [500000][50]f32 -g [100][50]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk30d50.txt
	#futhark dataset -b -g 1i32 -g 60i64 -g [500000][50]f32 -g [100][50]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk60d50.txt
	#futhark dataset -b -g 1i32 -g 120i64 -g [500000][50]f32 -g [100][50]f32 > tmpdata
	#./$(EXEAcc) -t /dev/stdout < tmpdata >> super_acc_t1resk120d50.txt
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][150]f32 -g [100][150]f32 > tmpdata
	./$(EXEAcc) -t /dev/stdout < tmpdata >> acc_t1resk1d150.txt
	futhark dataset -b -g 1i32 -g 5i64 -g [500000][150]f32 -g [100][150]f32 > tmpdata
	./$(EXEAcc) -t /dev/stdout < tmpdata >> acc_t1resk5d150.txt
	futhark dataset -b -g 1i32 -g 15i64 -g [500000][150]f32 -g [100][150]f32 > tmpdata
	./$(EXEAcc) -t /dev/stdout < tmpdata >> acc_t1resk15d150.txt
	futhark dataset -b -g 1i32 -g 30i64 -g [500000][150]f32 -g [100][150]f32 > tmpdata
	./$(EXEAcc) -t /dev/stdout < tmpdata >> acc_t1resk30d150.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][150]f32 -g [100][150]f32 > tmpdata
	./$(EXEAcc) -t /dev/stdout < tmpdata >> acc_t1resk60d150.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][150]f32 -g [100][150]f32 > tmpdata
	./$(EXEAcc) -t /dev/stdout  < tmpdata >> acc_t1resk120d150.txt

compDC:
	futhark cuda knn-driver.fut
	futhark cuda knn-comparer.fut

runD: $(EXEDriver) $(EXEAcc)
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3 < tmpdata  >> super_t1resk1d15.txt
	futhark dataset -b -g 1i32 -g 5i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3 < tmpdata  >> super_t1resk5d15.txt
	#./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk5d15.txt
	futhark dataset -b -g 1i32 -g 15i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3 < tmpdata >> super_tt1resk15d15.txt
	#./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk15d15.txt
	futhark dataset -b -g 1i32 -g 30i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n  -r 3 < tmpdata >> super_tt1resk30d15.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3 < tmpdata  >> super_t1resk60d15.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3 < tmpdata  >> super_t1resk120d15.txt
	#./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk30d15.txt
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_tt1resk1d50.txt
	futhark dataset -b -g 1i32 -g 5i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_tt1resk5d50.txt
	#./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk5d50.txt
	futhark dataset -b -g 1i32 -g 15i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_tt1resk15d50.txt
	#./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk15d50.txt
	#./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t1resk30d50.txt
	futhark dataset -b -g 1i32 -g 1i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_t1resk1d150.txt
	futhark dataset -b -g 1i32 -g 5i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_t1resk5d150.txt
	futhark dataset -b -g 1i32 -g 30i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_tt1resk30d50.txt
	futhark dataset -b -g 1i32 -g 15i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_t1resk15d150.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_tt1resk60d50.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_tt1resk120d50.txt
	futhark dataset -b -g 1i32 -g 30i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_t1resk30d150.txt
	futhark dataset -b -g 1i32 -g 60i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_t1resk60d150.txt
	futhark dataset -b -g 1i32 -g 120i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 3  < tmpdata >> super_t1resk120d150.txt


runDt10: $(EXEDriver)
	futhark dataset -b -g 10i32 -g 5i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2 < tmpdata >> t10resk5d15.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk5d15.txt
	futhark dataset -b -g 10i32 -g 15i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk15d15.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk15d15.txt
	futhark dataset -b -g 10i32 -g 30i64 -g [500000][15]f32 -g [500000][15]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2 < tmpdata >> t10resk30d15.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk30d15.txt
	futhark dataset -b -g 10i32 -g 5i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk5d50.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk5d50.txt
	futhark dataset -b -g 10i32 -g 15i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk15d50.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk15d15.txt
	futhark dataset -b -g 10i32 -g 30i64 -g [500000][50]f32 -g [500000][50]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk30d50.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk30d50.txt
	futhark dataset -b -g 10i32 -g 5i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk5d150.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk5d150.txt
	futhark dataset -b -g 10i32 -g 15i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk15d150.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk15d150.txt
	futhark dataset -b -g 10i32 -g 30i64 -g [500000][150]f32 -g [500000][150]f32 > tmpdata
	./$(EXEDriver) -t /dev/stdout -n -r 2  < tmpdata >> t10resk30d150.txt
	./$(EXEAcc) -t /dev/stdout < tmpdata  >> acc_t10resk30d150.txt


runtests:
	make runD
	make runDt10
	make runB

clean: 
	rm -f *.c
	rm -f *.actual
	rm -f *.expected
