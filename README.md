# AntVenture - README
Isaac Planas-Sitjà

2024-06-04

## AntVenture description
This is a description of the software AntVenture, which was created for education and research purposes. It requires PySimpleGUI, so you will be asked to create an account. This software simulates a colony of ants foraging on a liquid food source. AntVenture has been coded using real experimental data in Diacama ants (see Fujioka etal paper).

Using AntVenture, you can investigate the optimal foraging strategy depending on the strategy used by ants (trophallaxis vs social bucket) in diverse social and environmental contexts. You can modify the number of foraging ants and size of the colony, as well as the distance between the colony and the food source, difficulty of the terrain and quality of the food (= quantity of sugar; yeah, ants are sweet, they love sugar!). With AntVenture you will measure how the food is spread inside the colony and how fast ants are fed depending on the parameters you choose. **The data that you will obtain will be very similar to what you would obtain if you were conducting real experiments with ants!** This program is intended for students to: a) practice statistics, b) practical tool to learn more about ants social behaviour and concepts related to social feedback or sociality, c) simulate real data when scholar do not have the time or the means to conduct experiments with ants, d) serve as a basic tool with customizable parameters to potentially investigate other social behaviours in ants or other social insects.

### Trophallaxis vs Social Bucket
Ants mainly possess two strategies to collect food in liquid form. They can either “drink” the liquid and store it in their crop, or they can “grab” a drop of liquid with their mandibules, and transport it as such (or they can use both strategies at the same time!). When drinking, ants will use what is called **trophallaxis** to share the liquid with other ants in the colony, where foragers regurgitate this fluid to distribute to nest-mates. When using grabbing, ants will use what is called the **social bucket** to share food: other ants directly drink from the drop between their mandibles when they are inside the colony. When using both strategies, well… first they'll need to drink and grab the liquid food, go to the nest, and share it first using the social bucket strategy, and then start doing trophallaxis.

| ![Ken Naganawa image](https://github.com/iplanass/AntVenture/blob/main/images/im1.png "Ken Naganawa image") |
|:--:|
| Image from Fujioka, H., Marchand, M. & LeBoeuf, A. C. Diacamma ants adjust liquid foraging strategies in response to biophysical constraints. Proc. R. Soc. B Biol. Sci. 290, (2023). |

While both strategies are used, one strategy may be better than the other depending on the context. For instance, grabbing a drop of liquid is very fast, independently of the viscosity of the liquid, but there is a chance of dropping the liquid on the floor when returning to the nest. On the contrary, the speed of *drinking* decreases with the increase of viscosity, but it is the safest way to transport the liquid back to the nest. Once forager ants arrive to the nest, they will pass this liquid to other ants, and will stimulate other workers to go out and forage on that food source. The speed to pass the liquid to other ants will highly depend on the strategy used, trophallaxis, social bucket, or both.

### Select working directory

![alt text](https://github.com/iplanass/AntVenture/blob/main/images/im2.png "working directory")

Introduce your working directory where you want to save the results. The name of the results will be automatically generated using date and time (this is to avoid replacing previous results using the same name by mistake).

## Parameters

![alt text](https://github.com/iplanass/AntVenture/blob/main/images/im3.png "parameter panel")

**Number of simulations:** How many simulations you want to run.

**Number of ants:** Size of the colony, counting foragers (extra-nidal ants, they go foraging) and nurses (intra-nidal ants, they stay inside the nest).

**Number of foragers:** This parameter defines the maximum numbers that can go out of the nest (nurses = Num of ants - Num foragers).

**Distance (food source):** Defines the distance between the colony and food source. Longer the distance, higher the chances to drop the liquid if using “grabbing”.

**Simulation time (s):** Time of the experiment in seconds.

![alt text](https://github.com/iplanass/AntVenture/blob/main/images/im4.png "simulation options")

**Sugar or Viscosity and value:** Choose whether the value you introduce is about Sugar concentration [0 - 1] or measured viscosity (in mPa; the maximum value is ~0.054).

**Terrain difficulty:** Choose the difficulty of the terrain (by default = 0). Zero means that the terrain is flat and smooth, so ants have no difficulty in walking between the nest and the food source. Increasing the terrain difficulty will increase the difficulty for ants using “social bucket” (or “both” strategy), as they will have higher chances of loosing the droplet. This parameter has no effect if ants are using trophallaxis.

![alt text](https://github.com/iplanass/AntVenture/blob/main/images/im5.png "select strategy")

**Method:** this method affects only the social bucket and both strategies. Simple method uses a simpler (thus faster) way to compute how food is spread inside a colony compared to complex. For details see below.

Finally, choose an option, whether your ants are using Trophallaxis, Social Bucket of Both strategies. Remember to click SAVE if you want to save the parameters (otherwise they will be empty next time you open the app). Then, just click RUN, and have fun!

### Simple vs Complex:

When ants use social bucket, they can feed several ants at a time (~1 - 4 nest-mates) when they are in the nest. There is a chance that some of these ants are overfed, that means they take more liquid than they need (this has been empirically observed). In such case, they can pass this remaining liquid to other hungry ants in the nest. When we allow them to do that, is what we call the complex method. When choosing simple method, we kind of… simplify and speed up things. In other words, ants that are overfed do not continue to pass this liquid, we just consider that they had an overdose of sugar. While they use slightly different mechanisms, given the stochasticity of the events, simulations with both methods will lead to very similar results at the end.

## Examples
Simulations using trophallaxis:

![alt text](https://github.com/iplanass/AntVenture/blob/main/images/im6.png "example trophallaxis")

Simulations using social bucket:

![alt text](https://github.com/iplanass/AntVenture/blob/main/images/im7.png "example social bucket")

### Notes
This software was tested on Windows 10, LinuxOS (Manjaro Xfce 4.20 kernel version 6.12) and MacOS 15.3.1 (it may not run in older MacOS versions).
