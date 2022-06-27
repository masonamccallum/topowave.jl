using Javis

function object(p=O, color="black")
    sethue(color)
    circle(p, 25, :fill)
    return p
end

begin
myvideo = Video(500, 500);
Background(1:70, ground);
red_ball = Object(1:70, (args...) -> object(O, "red"), Point(100,0));
blue_ball = Object(1:70, (args...) -> object(O, "blue"), Point(200,80));
act!(red_ball,Action(anim_rotate_around(2π, O)));
act!(blue_ball,Action(anim_rotate_around(2π, O.0,red_ball)));

render(
    myvideo;
    pathname="vis_play.gif"
)
end