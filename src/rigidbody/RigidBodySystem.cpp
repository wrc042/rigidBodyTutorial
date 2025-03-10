#include "rigidbody/RigidBodySystem.h"

#include "collision/CollisionDetect.h"
#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"

#include "solvers/SolverBoxPGS.h"

namespace Eigen
{
    typedef Matrix<float, 6, 1, ColMajor> Vector6f;
    typedef Matrix<float, 4, 3, ColMajor> KinematicMap;
}

namespace
{
    static inline Eigen::Quaternionf kinematicMap(const Eigen::Quaternionf& q, const Eigen::Vector3f& omega)
    {
        return 0.5f * Eigen::Quaternionf(
            -q.x() * omega.x() - q.y() * omega.y() - q.z() * omega.z(),
            q.w() * omega.x() + q.z() * omega.y() - q.y() * omega.z(),
            -q.z() * omega.x() + q.w() * omega.y() + q.x() * omega.z(),
            -q.y() * omega.x() - q.x() * omega.y() + q.w() * omega.z()
        );
    }
}

RigidBodySystem::RigidBodySystem() :
    contactStiffness(1e6f), contactDamping(1e5f),
    mu(0.4f), solverIter(10), m_preStepFunc(nullptr), m_resetFunc(nullptr)
{
    m_collisionDetect = std::make_unique<CollisionDetect>(this);
    m_solver = new SolverBoxPGS(this);
}

RigidBodySystem::~RigidBodySystem()
{
    clear();
}

void RigidBodySystem::addBody(RigidBody *_b)
{
    // Add the body.
    m_bodies.push_back(_b);
}

void RigidBodySystem::step(float dt)
{
    // Initialize the system.
    // Apply gravitional forces and reset angular forces.
    // Cleanup contacts from the previous time step.
    //
    for(auto b : m_bodies)
    {
        b->f = b->mass * Eigen::Vector3f(0.f, -9.81f, 0.f);
        b->tau.setZero();
        b->fc.setZero();
        b->tauc.setZero();
        b->contacts.clear();
    }

    // Standard simulation pipeline.
    computeInertias();

    if( m_preStepFunc )
    {
        m_preStepFunc(m_bodies);
    }

    // Detect collisions and compute the contact Jacobians.
    // The list of contacts is automatically cleared and regenerated by detectCollisions(),
    //  and computeContactJacobians() will traverse this list to build contact Jacobians.
    //
    m_collisionDetect->detectCollisions();
    m_collisionDetect->computeContactJacobians();

    // Update contact stiffness and damping for all contacts.
    //
    auto contacts = m_collisionDetect->getContacts();
    for(auto c : contacts)
    {
        c->k = contactStiffness;
        c->b = contactDamping;
        c->mu = mu;
    }

    for(auto b : m_bodies)
    {
        b->fc.setZero();
        b->tauc.setZero();
    }

    // Compute constraint forces.  
    // The LCP solver will be called here.
    //
    calcConstraintForces(dt);

    // TODORBT Perform numerical integration to first update the velocities of each rigid body in @a m_bodies, 
    //      followed by the positions and orientations.
    //
    for(RigidBody* b : m_bodies)
    {
        if (b->fixed) {
            b->xdot = Eigen::Vector3f::Zero();
            b->omega = Eigen::Vector3f::Zero();
            continue;
        }
        // velocity
        b->xdot += dt * (1.0f / b->mass) * (b->f + b->fc);
        b->omega += dt * b->Iinv * (b->tau + b->tauc - b->omega.cross(b->I * b->omega));

        // position
        b->x += dt * b->xdot;

        // rotation
        b->q.coeffs() += dt * kinematicMap(b->q, b->omega).coeffs();
        b->q.normalize();
    }
}

void RigidBodySystem::clear()
{
    if( m_resetFunc )
    {
        m_resetFunc();
    }

    m_collisionDetect->clear();

    for(auto b : m_bodies)
    {
        delete b;
    }
    m_bodies.clear();
}



void RigidBodySystem::computeInertias()
{
    for(RigidBody* b : m_bodies)
    {
        b->updateInertiaMatrix();
    }
}

// Accessors for the contact constraint array.
const std::vector<Contact*>& RigidBodySystem::getContacts() const
{
    return m_collisionDetect->getContacts();
}

std::vector<Contact*>& RigidBodySystem::getContacts()
{
    return m_collisionDetect->getContacts();
}

void RigidBodySystem::calcConstraintForces(float dt)
{
    // Solve for the constraint forces lambda
    //
    m_solver->setMaxIter(solverIter);
    m_solver->solve(dt);

    // Apply the constraint forces as forces and torques acting on each body.
    // Essentially, we compute contact forces by computing :
    //
    //       f_contact = J^T * lambda
    //
    // for each contact constraint.
    //
    auto contacts = m_collisionDetect->getContacts();
    for(const auto c : contacts)
    {
        // Convert impulses in c->lambda to forces.
        //
        const Eigen::Vector6f f0 = c->J0.transpose() * c->lambda / dt;
        const Eigen::Vector6f f1 = c->J1.transpose() * c->lambda / dt;

        if (!c->body0->fixed)
        {
            c->body0->fc += f0.head<3>();
            c->body0->tauc += f0.tail<3>();
        }
        if (!c->body1->fixed)
        {
            c->body1->fc += f1.head<3>();
            c->body1->tauc += f1.tail<3>();
        }
    }
}

