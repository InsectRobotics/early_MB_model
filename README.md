# early_MB_model
Code relating to the original (2016)  MB model for route following using a spiking neuron model
These files were recovered from the SVN repository used by the insect robotics group prior to github. They were deposited on 10/10/2013 by Fei Peng, not updated since.

They relate to this paper: 
Ardin, Paul, Fei Peng, Michael Mangan, Konstantinos Lagogiannis, and Barbara Webb. "Using an insect mushroom body circuit to encode route memory in complex natural environments." PLoS computational biology 12, no. 2 (2016): e1004683.

1. Shang's model with a size of 49/360 (PN/KC) - the first model, done as a student project
2. Fei's CPU version (running on CPU), with a size of 360/4000(PN/KC), good
 for demonstrating and understanding the network behaviour since there are several recordings and plottings;
3. Fei's GPU version (running with MATLAB CUDA, asking for nvidia graphic 
card), which doesn't have plottings for network behaviour in a single run, but does have recordings and figures for overall results after training with a certain 
number of images. 
