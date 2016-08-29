
make 目录下共有3个文件夹：

1. ivf: 采用 Intel Visual Fortran 编译器，适用于Windows XP和Windows 7等操作系统；
2. cvf：采用 Compaq Visual Fortran 编译器，仅用于Windows XP操作系统；
3. gnu：采用GNU Fortran 编译器，用于linux、Mac OS操作系统；

注意事项：
1.  为了避免出现错误"error #6633: The type of the actual argument differs from the type of the dummy argument."，
    需关闭编译器的函数接口类型检查功能。对于Intel Visual Fortran编译器，可在属性页中将
    "Fortran|Diagnostics|Language Usage Warnings|Check Routine Interfaces"选项的值设为"No"；
2.  将属性页中的"Debugging|Working Directory"的值设为工作目录(即输入数据文件"stap90.in"所在的目录)，
    如设为"D:\TechPro\stap90\data"。程序将从该目录中读入stap90.in文件，并在该目录中输出结果数据文件。