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

# Taking it further: Checking the stability of the evolutionarily singular strategy

The script ```solvedbo.m``` above finds the evolutionarily singular strategy, however ideally we'd like to demonstrate that this singular strategy is both ESS-stable (meaning that no nearby variant strategy can invade it) and that it is convergent stable (meaning that it is an evolutionary attractor).

After running the script

    > solvedbo;

the workspace will now have access to the parameter dictionary ```p```, a vector of the optimal hatching times $x_c$

    x_cV =

    Columns 1 through 9:
    175.00   174.08   173.16   172.24   171.33   170.41   169.49   168.57   167.65
    ...
    etc.

and the corresponding singular $(y^*,z^*)$

    yzV =

    130.3250    30.6742
    130.3282    29.7526
    130.3281    28.8343
    130.3279    27.9160
    ...
    etc.

To check that each singular strategy is ESS-stable

    > check_ess(p,yzV,x_cV)

will return the eigenvalue of the Jacobian of the selection gradient with largest real part, which must be less than 0. To check that each singular strategy is convergent stable

    > check_conv(p,yzV,x_cV)

will return the eigenvalue of the Hessian with largest real part, which must be less than 0.

Refer to the Supplementary Material for explanations of the relationships between stability constraints and the eigenvalues.

