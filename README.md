# Reference

Kristensen, N.K., Johansson, J., Ripa, J. and Jonzen, N. Phenology of arrival and laying time in migratory birds in response to climate change. Paper in preparation.

# Quick Start

If you have git:

    $ git clone https://github.com/nadiahpk/phenology-two-trait-migratory-bird
    $ cd phenology-two-trait-migratory-bird

If you don't have git, just download and unpack the latest zip
file:

[https://github.com/nadiahpk/phenology-two-trait-migratory-bird/archive/master.zip](https://github.com/nadiahpk/phenology-two-trait-migratory-bird/archive/master.zip)

In Octave

    > solvedbo;

will produce the file ```u_q125.dat```. You can then run the gnuplot script

    $ gnuplot u_q125.gnu

to produce the figure file ```u_q125.eps```, which is the first pane in Figure 1.

![Figure 1a](https://raw2.github.com/nadiahpk/phenology-two-trait-migratory-bird/1d687030b0f4703d820244f8209a07b0d38f4c69/u_q125.png)

# Taking it further

## Changing parameter values

The function ```solvedbo``` can also be used to explore how Figure 1 changes with different parameter values. For example, the quickstart above could also be generated in the following way.

First, create a dictionary ```p``` with your default parameter values. I have stored one in the script ```params.m``` that you can use

    > params;
    > p
    p =

        scalar structure containing the fields:

            K =  100
            b_e =  0.030000
            a =  3
            sigma =  5
            Q_0 =  0.050000
            b_q =  0.50000
            u_q =  134
            L_m =  170
            b_s =  0.075000
            M =  0.50000
            z_f =  40
            z_n =  14

Let's say I'd like to look at the effect of changing the parameter ```u_q``` to 125. I want to plot the evolutionarily singular phenology over the optimal hatching time range of 135 - 175. I have a fair idea that the arrival date when optimal hatching is late is about day 120 (implying an approximate prelaying period of 175-120-z_n), and I'd like to print the results to a .dat file. I would run

    > [x_cV,yzV,nV]=solvedbo(p,'u_q',125,175,135,120,175-120-p.z_n,1)

This will write a .dat file ```u_q125.dat``` as before. I will also now have the optimal hatching times ```x_cV``` and evoluttionarily singular strategies ```yzV``` in the workspace

    x_cV =

    Columns 1 through 6:
    175.00   174.18   173.37   172.55   171.73   170.92
    ...
    etc.

## Checking the stability of the evolutionarily singular strategy

The script ```solvedbo.m``` above finds the evolutionarily singular strategy, however ideally we'd like to demonstrate that this singular strategy is both ESS-stable (meaning that no nearby variant strategy can invade it) and that it is convergent stable (meaning that it is an evolutionary attractor).

After running the function above, you will have ```p```, ```x_cV```, and ```yzV``` in the workspace. NOTE if you've used the ```solvedbo``` function the output will not match ```p``` for the change-variable, so update that first

    > p.u_q = 125

Now we can test evolutionary stability. There are two types of stability -- ESS stability and convergence stability -- and there are two sets of code to check for each -- one numerical and one analytic. Let's use check the ESS-stability of our system numerically

    > [numess,numHess] = numcheck_ess(p,yzV,x_cV);
    > numess
    numess =
    -3.8428e-04
    -3.8418e-04
    ...
    etc.

The vector ```numess``` is the eigenvalue of the Hessian matrix with the largest real part, as found numerically. ```numHess``` is the Hessian matrix of the last evolutionarily singular strategy evaluated, so you can use this to take a closer look at single points. In order to establish ESS-stability, this eigenvalue must be less than zero (see paper). You could plot this ```plot(x_cV,numess)``` to see how it changes with optimal hatching time.

Let's check the convergence stability of our system semi-analytically

    > [conv,Jac] = check_conv(p,yzV,x_cV);
    > conv
    conv =
    -3.8429e-04
    -3.8418e-04
    ...
    etc.

The vector ```conv``` is the is the eigenvalue of the Jacobian matrix with the largest real part, and it must also be less than zero in order to establish convergence stability. If you open ```check_conv.m``` you can see expressions for the derivatives of each of the components. These were found symbolically using Sage. The function ```numcheck_conv.m``` was used to verify numerically that the expressions were correct, and likewise there is a ```check_ess.m``` which provides the semi-analytic counterpart to ```numcheck_ess.m``` above.

Refer to the Supplementary Material of the paper for further explanations of the relationships between stability constraints and the eigenvalues.

# License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this software, either in source code form or as a compiled binary, for any purpose, commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this software dedicate any and all copyright interest in the software to the public domain. We make this dedication for the benefit of the public at large and to the detriment of our heirs and successors. We intend this dedication to be an overt act of relinquishment in perpetuity of all present and future rights to this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to
<http://unlicense.org/>
