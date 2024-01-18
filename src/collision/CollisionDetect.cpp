#include "collision/CollisionDetect.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

CollisionDetect::CollisionDetect(RigidBodySystem* rigidBodySystem) : m_rigidBodySystem(rigidBodySystem)
{

}

void CollisionDetect::detectCollisions()
{
    // First, clear any existing contacts.
    //
    clear();

    // Next, loop over all pairs of bodies and test for contacts.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for(unsigned int i = 0; i < bodies.size(); ++i)
    {
        for(unsigned int j = i+1; j < bodies.size(); ++j)
        {
            RigidBody* body0 = bodies[i];
            RigidBody* body1 = bodies[j];

            // Special case: skip tests for pairs of static bodies.
            //
            if (body0->fixed && body1->fixed) 
                continue;

            // Test for sphere-sphere collision.
            if( body0->geometry->getType() == kSphere &&
                body1->geometry->getType() == kSphere )
            {
                collisionDetectSphereSphere(body0, body1);
            }
            // Test for sphere-box collision
            else if( body0->geometry->getType() == kSphere &&
                     body1->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body0, body1);
            }
            // Test for box-sphere collision (order swap)
            else if( body1->geometry->getType() == kSphere &&
                     body0->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body1, body0);
            }
        }
    }
}

void CollisionDetect::computeContactJacobians()
{
    // Build constraint Jacobians for all contacts
    //
    for(auto c : m_contacts)
    {
        c->computeContactFrame();
        c->computeJacobian();
    }
}

void CollisionDetect::clear()
{
    // First, remove all contacts from rigid bodies.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for (auto b : bodies)
    {
        b->contacts.clear();
    }

    // Then, cleanup the local contact array.
    //
    for(auto c : m_contacts)
    {
        delete c;
    }
    m_contacts.clear();

}

void CollisionDetect::collisionDetectSphereSphere(RigidBody* body0, RigidBody* body1)
{
    Sphere* sphere0 = dynamic_cast<Sphere*>(body0->geometry.get());
    Sphere* sphere1 = dynamic_cast<Sphere*>(body1->geometry.get());

    // Implement sphere-sphere collision detection.
    // The function should check if a collision exists, and if it does
    // compute the contact normal, contact point, and penetration depth.
    //
    Eigen::Vector3f vec = body0->x - body1->x;

    const float rsum = (sphere0->radius + sphere1->radius);
    const float dist = vec.norm();
    if( dist < rsum )
    {
        const Eigen::Vector3f n = vec / dist;
        const Eigen::Vector3f p = 0.5f * ((body0->x - sphere0->radius*n) + (body1->x + sphere1->radius*n));
        const float phi = dist-rsum;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}

void CollisionDetect::collisionDetectSphereBox(RigidBody* body0, RigidBody* body1)
{
    // TODORBT Implement sphere-box collision detection.
    //      The function should check if a collision exists.
    // 
    //      If it does, compute the contact normal, contact point, and penetration depth and
    //      create a Contact and add it to m_contacts.
    //

    Sphere* sphere = dynamic_cast<Sphere*>(body0->geometry.get());
    Box* box = dynamic_cast<Box*>(body1->geometry.get());

    Eigen::Vector3f pos_loc = body1->q.inverse() * (body0->x - body1->x);

    Eigen::Vector3f g;
    Eigen::Vector3f bdim_half = box->dim * 0.5f;
    float dphi = std::numeric_limits<float>::max();
    for(int i=0; i < 3; i++){
        float v = pos_loc[i];
        v = (std::max)(v, -bdim_half[i]);
        v = (std::min)(v, bdim_half[i]);
        g[i] = v;
    }

    if((g - pos_loc).norm() > 1e-6){
        Eigen::Vector3f normal = Eigen::Vector3f::Zero();
        float phi = (g - pos_loc).norm() - sphere->radius;
        if(phi < 0){
            normal = body1->q * (pos_loc - g);
            normal.normalize();
            m_contacts.push_back( new Contact(body0, body1, body1->q * g + body1->x, normal, phi) );
        }
    }else{
        float phi = std::numeric_limits<float>::max();
        Eigen::Vector3f normal = Eigen::Vector3f::Zero();
        for(int i = 0; i < 3; i++){
            if(abs(bdim_half[i] - pos_loc[i]) < phi){
                phi = abs(bdim_half[i] - pos_loc[i]);
                normal = Eigen::Vector3f::Unit(i);
                g = pos_loc + normal * phi;
            }
            if(abs(bdim_half[i] + pos_loc[i]) < phi){
                phi = abs(bdim_half[i] + pos_loc[i]);
                normal = -Eigen::Vector3f::Unit(i);
                g = pos_loc + normal * phi;
            }
        }
        phi = -phi;
        m_contacts.push_back( new Contact(body0, body1, body1->q * g + body1->x, normal, phi) );
    }
}
