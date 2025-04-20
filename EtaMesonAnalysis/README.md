This section (folder) contains code used for Eta Meson Analysis. We also use it to reconstruct the background via event mixing. As of Sept 9, 2024, we are still working on to get background via position / momentum swapping and have not yet updated here.


Update:: April, 20 2024

There are two versions of event mixing. The one with "OLD" suffix is one that was in existence since begining and is computationally inefficient process. The other one is much better (produce high stat backgroun) and computationally more efficient process. Both function have same algorithm, thus make same signal (foreground as expected) and same nature of background (but different statistics). We update the trigger as well (in all functions) and also in process event. There is upate in macro. We do not need TChain macro in addition and all input parameters in .cc file will be passed from macro thus removing the requirement of building for every parameter change.

As of today, we are still working on position-swapping method so the current version in code may not work. Once updated/tested, we will upload in github and update this README file.
