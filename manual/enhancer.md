# Evaluation of transgenic fish for enhancer activity
*Methods for processing images and analyzing GFP intensity scores from transgenic fish*

**Relevant Figure**: Figure 4

## Adjust color and brighness in images
Apply identical brightness and contrast shifts to images shown in Figure 4B
```
sbatch scripts/enhancer/adjust_brightness_contrast.sh
```

## Plot GFP intensity
Plot normalized GFP intensity scores and perform statistical test:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/enhancer_figure.Rmd](../scripts/R/enhancer_figure.Rmd)  
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/enhancer_figure.html](../scripts/R/enhancer_figure.html)  
