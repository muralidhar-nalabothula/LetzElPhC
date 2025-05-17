#include "print_letz_logo.h"

#include <stdio.h>

void print_letz_logo(FILE* output)
{
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
}
