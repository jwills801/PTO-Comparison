# READ ME


OSWEC_Hydraulic_PTO.slx - This code used to simulate the whole check valve PTO. There is an example file from wecSim PTOSIM which does this. I had
adjusted this simulink model so that it uses our pump map for the motor and our 90% efficiency for the generator. It was hard to set parameters
for this case - you had to make sure that the volume of the accumulators ended the simulation at the same value it began with. Also how large should the accums
be? What should the precharge pressure be?

For these reasons, and to match what we did for the HHEA case, the accumulators were assumed to be large enough to make any volume deviation negligible. 
There was a constant pressure difference Between the two accums. This pressure difference is optimized via a grid search. The volume constaint is met by
correctly sizing the motor.