#ifndef __GNUPLOT_HPP
#define __GNUPLOT_HPP

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>
#include <fcntl.h>

/** Maximal number of simultaneous temporary files */
#define GP_MAX_TMP_FILES    64
/** Maximal size of a temporary file name */
#define GP_TMP_NAME_SIZE    512

/*---------------------------------------------------------------------------
                                New Types
 ---------------------------------------------------------------------------*/

typedef struct _GNUPLOT_CTRL_ {

  /** Pipe to gnuplot process */
  FILE    * gnucmd ;

  /** Number of currently active plots */
  int       nplots ;

  /** Current plotting style */
  char      pstyle[32] ;

  /** Name of temporary files */
  char      to_delete[GP_MAX_TMP_FILES][GP_TMP_NAME_SIZE] ;

  /** Number of temporary files */
  int       ntmp ;

} gnuplot_ctrl;

/* Class declaration */
class GnuPlot {

protected:

private:

  gnuplot_ctrl * plot_ctrl;

public:

  GnuPlot();
  ~GnuPlot();

  char * gnuplot_get_program_path(char * pname);
  gnuplot_ctrl * gnuplot_init(void);
  void gnuplot_close(void);

  void gnuplot_cmd( char *  cmd, ...);
  void gnuplot_setstyle(char * plot_style);
  void gnuplot_set_xlabel(char * label);
  void gnuplot_set_ylabel(char * label);
  void gnuplot_set_title(char * title);
  void gnuplot_resetplot(void);
  void gnuplot_plot_x(double * d, int n, char * title);
  void gnuplot_plot_xy( double *x, double *y, int n, char *title );
  void gnuplot_plot_once( char *title, char *style, char *label_x,
			 char *label_y, double *x, double *y, int n );
  void gnuplot_plot_slope( double a, double b, char *title );
  void gnuplot_plot_equation(char * equation, char * title) ;
  void gnuplot_plot_histogram( double *ordinate, double *rawdata, 
			       const unsigned int nbins,
			       int overflow, char * title
			       );
};
#endif
#endif
