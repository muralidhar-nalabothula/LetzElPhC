/**
 * @file
 * @brief Implementation of functions for printing informational messages and
 * logos.
 *
 * This file contains the implementation of output functions, primarily used
 * for displaying the program's logo and version information to the user.
 */
#include <stdio.h>

#include "print_info.h"

/**
 * @brief Prints the ELPH logo to the specified output stream.
 *
 * The logo is printed only by the MPI process with rank 0 to avoid
 * duplicate output in parallel runs.
 *
 * @param mpi_rank The MPI rank of the calling process. The logo is printed if
 * this is 0.
 * @param output A pointer to the @c FILE stream (e.g., @c stdout or a log file)
 * where the logo should be written.
 */
void print_ELPH_logo(int mpi_rank, FILE* output)
{
    if (mpi_rank)
    {
        return;
    }
    fprintf(output,
            "       */*                                                        "
            "    *+        \n");
    fprintf(output,
            "      */                                                          "
            "     //       \n");
    fprintf(output,
            "    +/*   +*+            +          +++++++ ++   +++++   ++      "
            "+/***/  */+    \n");
    fprintf(output,
            "   */+    +*+     +++++ ***+ ++++*  ++      ++   ++   ++ ++++++  "
            "/*   +    /*   \n");
    fprintf(output,
            "  //      +*+    ** ++*  *+    +*   +++++++ ++   ++++++  ++  ++ "
            "+/+         */+ \n");
    fprintf(output,
            "  */+     +*+    **      *+   *+    ++      ++   ++      ++  ++  "
            "/*   *+    */* \n");
    fprintf(output,
            "   +/*     +++++  +++++  +++ +++++  +++++++ ++   ++      ++  ++   "
            "+***+    */   \n");
    fprintf(output,
            "     /*                                                           "
            "      +/*     \n");
    fprintf(output,
            "      */+            **//*+        **+     **        +**/*+       "
            "     */+      \n");
    fprintf(output,
            "       +/+         +//****//      /  /+   /  /      ////////+     "
            "    //        \n");
    fprintf(output,
            "                   +//+++*//+++++     /  /+  +*+*+++////////+     "
            "              \n");
    fprintf(output,
            "                    +/////*+          +**+           */////+      "
            "              \n");
    fprintf(output, "\n\n");
}
