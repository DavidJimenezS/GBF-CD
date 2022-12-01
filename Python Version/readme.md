# GBF-CD Python Version
Graph based fusion for change detection of the "[Graph-Based Data Fusion Applied to: Change Detection and Biomass Estimation in Rice Crops](https://www.mdpi.com/2072-4292/12/17/2683)" published in MDPI.

Here you will find all the code used to generate a change map in fourteen cases of study.

generate_nystrom_graph.py is the main code that generates the changes map, and 
you can find the implementation in the main.py script.

Please if you use the datasets and/or the code cite us as:<br/>

@article{&nbsp;&nbsp;&nbsp;JimenezSierra2020graph,<br/>
         &nbsp;&nbsp;&nbsp;title={Graph-Based Data Fusion Applied to: Change Detection and Biomass Estimation in Rice Crops},<br/>
         &nbsp;&nbsp;&nbsp;author={Jimenez-Sierra, David Alejandro and Ben{\\'i}tez-Restrepo, Hern{\\'a}n Dar{\\'i}o and Vargas-Cardona, Hern{\\'a}n Dar{\\'i}o and Chanussot, Jocelyn},<br/>
         &nbsp;&nbsp;&nbsp;journal={Remote Sensing},<br/>
         &nbsp;&nbsp;&nbsp;volume={12},<br/>
         &nbsp;&nbsp;&nbsp;number={17},<br/>
         &nbsp;&nbsp;&nbsp;pages={2683},<br/>
         &nbsp;&nbsp;&nbsp;year={2020},<br/>
         &nbsp;&nbsp;&nbsp;publisher={Multidisciplinary Digital Publishing Institute}<br/>
        }
        
In addition, you can check the [Requirement](https://github.com/DavidJimenezS/GBF-CD/blob/master/Python%20Version/requirements.txt)) file to create an enviroment that runs the code as follows:

"conda create --name <env> --file requirements.txt"
