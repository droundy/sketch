# sketch-space-time

This is a drawing and animation program for children (specifically for one five-year-old child),
with an emphasis on conceptual ease.

## Getting started

You can run the program using `cargo run --release`, and can optionally specify a save file name
with `cargo run --release FILENAME`.  You can exit the program by hitting the escape key, which is
the only keyboard command.  The drawing is saved on exit to an internal (and not necessarily stable)
format.

Tools are selected in the upper left.

The color selector should be obvious.  It doesn't allow you
to select any possible color, but is pretty easy to use, and includes all the colors my daughter
is likely to want (and is easier to use than a full three dimensional selector).

Layers can be created, selected, and reordered in the lower left.  Each layer has a pen color and a fill
color (initially the fill is transparent) that can be modified.

The animation control is on top.  You can create new keyframes, move keyframes on the timeline, and start
or stop the animation itself.

## UX Principles

Actions should be as concrete as possible, and each action should be reversible by direct action.
e.g. if something is drawn, this should be reversible by erasing it.  If something is erased, undoing
that action should just require drawing it again.  This principle requires *layers*, since otherwise
if you drew on top of something, the reverse wouldn't be to erase that, but instead to erase it *and
then draw whatever you overwrote*.  In this program, each color (of a finite set of colors) is a
layer, so creating a new color always requires creating a new layer.

## UX Problems

The animation interface can be tricky to use for a 5-year old.

A keyframe once created cannot be deleted, which violates the principle of actions being reversible.
I'm not sure of a good UI for deletion, and deletion of course tends to be hard to reverse.  If there
become too many keyframes, the whole business becomes hard to manage.

The layers have a similar challenge, in that once created a layer cannot be deleted.  However, it is
relatively easy to change the colors and redefine a layer.

### Naming

This crate could really use a better name.  I will appreciate any suggestions!