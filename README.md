# MIC

Matlab Inversion Collaboration, including iterative methods and global optimization methods.

## License

[The MIT License](https://tchel.mit-license.org)

## Inversion methods

### Iterative methods
- [KA](Kaczmarz.m): Kaczmarz's Algorithm
- [ART](ART.m): Algebraic Reconstruction Technique
- [SIRT](SIRT.m): Simultaneous Iterative Reconstruction Technique
- [CG](ConjugateGradient.m): Conjugate Gradient method
- [LM](LevenbergMarquardt.m): Levenberg-Marquardt method
### Global optimization methods

- [GS](HighorderGridsearch.m): Grid-Search method
- [MC](MentoCarlo.m): Monte-Carlo method
- SA: Simulated Annealing method, including [Metropolis](Metropolis.m) and [HeatBath](HeatBath.m)
- [GA](StandardGeneticAlgorithm.m): Genetic Algorithm

## Example

There are two examples to demonstrate how to use these matlab function:

- [example1.m](example1.m) for LM, GS, MC, SA & GA
- [example2.m](example2.m) for KA, ART, SIRT & CG

## Reference

- Aster R C, Borchers B, Thurber C H. **Parameter estimation and inverse problems** [M]. Academic Press, 2011.
- Sen M K, Stoffa P L. **Global optimization methods in geophysical inversion** [M]. Cambridge University Press, 2013.

## Author

Tche LIU, seistche@gmail.com, USTC
