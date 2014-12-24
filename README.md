###ON CPU

>####require
+ CMake version 2.6 or newer
+ boost 1.55.0
+ g++ 4.7 or newer

> ####setup
~~~
$cmake .
$make
~~~
> ####usage
> build index
~~~
$./sbwt build -i [reference.fa] -p [ouput prefix]
~~~
> for more detail please use :
~~~
$./sbwt build -h
~~~

> alignment
~~~
$./sbwt map -i [reads.fq] -p [index prefix] -o [ouput.sam]
~~~
> for mor detail please use :
~~~
$./sbwt map -h
~~~
