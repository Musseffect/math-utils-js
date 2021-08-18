import { Epsilon, near } from "./utils";
import vec3 from "./vec3";

export default class axisAngle{
    axis:vec3;
    angle:number;
    constructor(axis:vec3, angle:number){
        this.axis = axis.normalized();
        this.angle = angle;
    }
    static near(a:axisAngle, b:axisAngle, threshold?:number):boolean{
        if(!threshold)
            threshold = Epsilon;
        return vec3.near(a.axis, b.axis, threshold) && near(a.angle, b.angle, threshold);
    }
    rotate(point:vec3):vec3{
        return vec3.add(vec3.lerp(this.axis.scale(vec3.dot(this.axis, point)), point, Math.cos(this.angle)),
         vec3.cross(this.axis, point).scaleSelf(Math.sin(this.angle)));
    }
    toString():string{
        return `{ axis:${this.axis.toString()}, angle:${this.angle.toFixed(4)} }`;
    }
    // TODO: euler angles to axis angle
}