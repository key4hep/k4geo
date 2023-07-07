CLD
====

FCCee_o2_v01
------------

This option was used for flavourtagging studies with the smaller radius beampipe and vertex detector, presented in the
(CLD Note)[https://arxiv.org/abs/1911.12230]. It reduces the radius and width of the two inner most barrel layers, and
preserves the length of the staves.

FCCee_o2_v02
------------

This model is an update of v01 with a fixed TrackerEndcapSupport

FCCee_o2_v03
------------

This model was taken from FCCDetectors `FCCee_o2_v02`. It doesn't fit in the option 2 category based on its too large
vertex detector radius. This model is _obsolete_.

FCCee_o2_v04
------------

This model was taken from FCCDetectors `FCCee_o2_v03`. It presents an alternate implementation of the vertex detector
with a smaller radius, which preserves angular coverage of the vertex barrel. Various issues of the Tracker structure
where inherited from the FCCDetectors `FCCee_o2_v02` model.

CLD_o2_v05
----------

This model is based on `FCCee_o2_v02` from k4geo. It has an updated design of the beampipe corresponding to the latest
standard design (low-impedance, no HOM absorber); changed Vertex Detector Barrel layers to fit the beampipe contraints (reduced barrel length); fixed
Overlaps in the Inner and Outer Tracker.


CLD_o3_v01
----------

CLD option 3 includes ARC detector for PID
