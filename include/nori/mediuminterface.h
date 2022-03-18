//
// Created by 郭彬 on 2022/3/18.
//

#ifndef NORI_MEDIUMINTERFACE_H
#define NORI_MEDIUMINTERFACE_H

#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

class Medium;

struct MediumInterface {
    MediumInterface(const Medium *medium)
            : m_inside(medium), m_outside(medium) {}

    MediumInterface(const Medium *inside, const Medium *outside)
            : m_inside(inside), m_outside(outside) {}

    bool isMediumTransition() const { return m_inside != m_outside; }

    const Medium  *m_inside;
    const Medium  *m_outside;
};

NORI_NAMESPACE_END

#endif //NORI_MEDIUMINTERFACE_H
