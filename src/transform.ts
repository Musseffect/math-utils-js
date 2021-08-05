import mat4 from "./mat4";
import quat from "./quat";
import { assert, near } from "./utils";
import vec3 from "./vec3";

export default class transform{
    translation:vec3;
    rotation:quat;
    scale:vec3;
    constructor(translation:vec3, rotation:quat, scale:vec3){
        this.translation = translation;
        this.rotation = rotation;
        this.scale = scale;
    }
    hasUniformScale(): boolean {
        return near(this.scale.x, this.scale.y) && near(this.scale.x, this.scale.z) && near(this.scale.y, this.scale.z);
    }
    clone():transform{
        return new transform(this.translation.clone(), this.rotation.clone(), this.scale.clone());
    }
    static identity(): transform{
        return new transform(new vec3(0, 0, 0), quat.identity(), new vec3(1, 1, 1));
    }
    toAffine():mat4{
        return mat4.fromTRS(this.translation, this.rotation, this.scale);
    }
    transformPoint(point:vec3):vec3{
        return this.rotation.rotate(vec3.mul(point, this.scale)).addSelf(this.translation);
    }
    transformVector(vector:vec3):vec3{
        return this.rotation.rotate(vec3.mul(vector, this.scale));
    }
    transformInversePoint(point:vec3):vec3 {
        return this.rotation.inverse().rotate(vec3.sub(point, this.translation)).divSelf(this.scale);
    }
    transformInverseVector(vector:vec3):vec3 {
        return this.rotation.inverse().rotate(vector).divSelf(this.scale);
    }
    // todo: test
    // only possible with uniform scale
    inverse():transform {
        assert(this.hasUniformScale(), "can't inverse transform with non-uniform scale");
        let invRot = this.rotation.inverse();
        return new transform(
            invRot.rotate(vec3.div(this.translation.inverse(), this.scale)),
            invRot,
            new vec3(1.0 / this.scale.x, 1.0 / this.scale.y, 1.0 / this.scale.z)
        ); 
    }
    toString(): string{
        return `{ t:${this.translation.toString()}, r:${this.rotation.toString()}, s:${this.scale.toString()} }`;
    }
}