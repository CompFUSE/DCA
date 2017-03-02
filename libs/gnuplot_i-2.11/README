
gnuplot_i: a short documentation
------------------------------------


gnuplot_i (formerly known as gnuplot_pipes) is yet another attempt to
provide a programmer-friendly set of routines to use gnuplot as a data
displayer from a C program.

The procedure to follow to display graphics in a gnuplot
session is:

    1.  Open a new gnuplot session, referenced by a handle
    of type (pointer to) gnuplot_ctrl. This is done by
    calling gnuplot_init():

    gnuplot_ctrl * h ;
    h = gnuplot_init() ;

    h will be used as the handle to the gnuplot session
    which was just opened, by all further functions.


    2.  Send optionally display configuration orders. The
    following functions are available:

    gnuplot_setstyle(handle, style)
        sets the plotting style of the next plots

    gnuplot_set_xlabel(handle, label)
        sets the X label for the next plots

    gnuplot_set_ylabel(handle, label) 
        sets the Y label for the next plots

    Examples:

    gnuplot_setstyle(h, "impulses") ;
    gnuplot_set_xlabel(h, "my X label") ;
    gnuplot_set_ylabel(h, "my Y label") ;


    The most useful routine is probably gnuplot_cmd(), which allows
    you to send character strings to gnuplot as though they were typed
    in by a user. This routine works in a printf fashion, accepting
    the same kind of format string and variable number of arguments.
    The prototype is:

    gnuplot_cmd(handle, command, ...)

    Example:

    char myfile[] = "/data/file_in.dat" ;
    int  i ;
    
    gnuplot_cmd(handle, "plot '%s'", myfile);
    for (i=0 ; i<10 ; i++) {
        gnuplot_cmd (handle, "plot y=%d*x", i);
    }


    The following commands request the output to be done to
    a postscript file named 'curve.ps':

    gnuplot_cmd(h, "set terminal postscript") ;
    gnuplot_cmd(h, "set output \"curve.ps\"") ;

    Using gnuplot_cmd(), it should be fairly easy to add up some
    more configuration-related functions whenever needed.


    3.  Send some display orders: either functions or
    doubles, or double points. The following functions are
    available:

    gnuplot_plot_slope()
        to display a slope

    gnuplot_plot_equation()
        to display an equation

    gnuplot_plot_x()
        to display user-defined 1d data with a single variable.
        Input is a list of values, assumed regularly spaced on the
        X-axis.

    gnuplot_plot_xy()
        to display user-defined 1d data with two variables.
        Input is two lists of point coordinates.

    gnuplot_resetplot()
        states that the current gnuplot display will be cleared
        before next display is done.



    4.  Close the gnuplot handle. This is very important,
    since otherwise all temporary files will NOT be cleaned
    for /tmp and /var/tmp.

    gnuplot_close(h) ;


See examples of gnuplot_i use in the provided files
'example.c' and 'anim.c'.


Some more points before you start using gnuplot_i
-------------------------------------------------



    gnuplot_i is completely free software. Use it for
	whatever you want to do with it without any fee, and do not
	hesitate to send feedback to me if you wish:
    
        <ndevilla@free.fr>

    If you can do it, I would appreciate a mention somewhere
    that you are using 'gnuplot_i' in your application.
    Something like:
    "This software uses the gnuplot_i library written by
    N.Devillard <ndevilla@free.fr>
    If you are using gnuplot_i for a web-based
    application, you can also add a link to the gnuplot home
    page:

        http://ndevilla.free.fr/gnuplot/

    If you are so happy about this piece of code that you would like
	to fund more stuff like that, do not hesitate to send me cash :-)


    Notice that it is possible to open several gnuplot sessions
    from the same program. Just be careful then to which session
    you are talking when using functions. Example:

    h1 = gnuplot_init() ;
    h2 = gnuplot_init() ;

    gnuplot_plot_equation(h1, "sin(x)", "sine on first session");
    gnuplot_plot_equation(h2, "log(x)", "log on second session") ;
    sleep(3) ;
    gnuplot_close(h1) ;
    gnuplot_close(h2) ;

    Never forget to close opened sessions! This would leave
    some temporary files. Hopefully, your temporary
    directories should be cleaned automatically every now
    and then, but best is simply not to forget to close all
    windows before leaving the house.

    User interactions are not part of gnuplot_i. Feel
    free to use your own.

    No timing mechanisms are provided to leave the outputs
    on screen until e.g. a key is pressed. I leave it up to
    gnuplot_i users to provide that kind of
    functionality in their own application, depending on the
    kind of interaction they need.

    No programming so far of 'splot' related commands.
    Should you implement some, I would recommend making use
    of the 'plot' related commands as a canvas.


To compile and run the examples
-------------------------------

You can use the provided Makefile, or compile by hand:

To compile the gnuplot_i module:
% cc -I. -c gnuplot_i.c

To compile both examples:
% cc -o example example.c gnuplot_i.o
% cc -o anim anim.c gnuplot_i.o

To compile if you want to use the Makefile:
% make
% make tests

Try it out to see if it works:
% test/anim				-- Will show an animated sine wave
% test/example			-- Will show various functions in use
% test/sinepng			-- Will generate a PNG file with a sine wave

The latter test requires gnuplot to have been compiled with PNG
support. This should be the case if you are using a pre-packaged
gnuplot under Linux (RPM).

Of course, you must have gnuplot installed and working on
your machine, accessible from the user account you are
using.

N. Devillard
